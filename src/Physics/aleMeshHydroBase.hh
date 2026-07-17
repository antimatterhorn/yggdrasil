// Copyright (C) 2026  Cody Raskin

#pragma once
#include "hydro.hh"
#include "../Mesh/alemesh.hh"
#include "../Boundaries/aleMeshBoundary.hh"
#include <iostream>
#include <algorithm>

// Forward declaration for return type of computeFlux (definition in HLL.cc,
// included by concrete solvers -- same pattern as GridHydroBase).
template<int dim>
struct HLLFlux;

// Finite-volume Eulerian hydro on Mesh::ALEMesh, walking Face objects instead
// of Grid's axis-aligned neighbor pairs. One NodeList entry per ALEMesh cell
// (cell-centered, positions taken from mesh->cellCentroid), same convention
// as GridHydroBase. Node positions are NOT moved by this class -- no
// mesh-velocity field, no remap; ALEMesh::updateGeometry() exists for a
// future pass that adds real node motion, not used here.
//
// Unlike GridHydroBase's EvaluateDerivatives (which divides by grid spacing,
// a shortcut that only works because Grid's face-area/cell-volume ratio is
// uniform everywhere), this does a real flux*area/volume divergence, since
// ALEMesh cells have no such uniformity to exploit.
template<int dim>
class ALEMeshHydroBase : public Hydro<dim> {
protected:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    Mesh::ALEMesh<dim>* mesh;
    mutable double dtmin = 1e30;

    // Which ALEMeshBoundary (if any) governs each face index, built once in
    // ZeroTimeInitialize from whichever registered boundaries are actually
    // ALEMeshBoundary instances.
    std::vector<ALEMeshBoundary<dim>*> faceBoundary;

public:
    ALEMeshHydroBase(NodeList* nodeList,
                      PhysicalConstants& constants,
                      EquationOfState* eos,
                      Mesh::ALEMesh<dim>* mesh)
        : Hydro<dim>(nodeList, constants, eos), mesh(mesh) {

        size_t numCells = mesh->getConnectivityMap().size();
        auto* pos = nodeList->getField<Vector>("position");
        for (size_t i = 0; i < numCells; ++i)
            pos->setValue(i, mesh->cellCentroid(i));

        State<dim>* state = &this->state;
        auto* v   = nodeList->getField<Vector>("velocity");
        auto* rho = nodeList->getField<double>("density");
        auto* u   = nodeList->getField<double>("specificInternalEnergy");

        state->template addField<Vector>(v);
        state->template addField<double>(rho);
        state->template addField<double>(u);
    }

    virtual ~ALEMeshHydroBase() {}

    virtual void
    ZeroTimeInitialize() override {
        EOSLookup();
        this->UpdateState();
        this->InitializeBoundaries();

        if (mesh->getFaces().empty())
            mesh->computeFaces();

        faceBoundary.assign(mesh->getFaces().size(), nullptr);
        for (auto* bc : this->boundaries) {
            auto* aleBc = dynamic_cast<ALEMeshBoundary<dim>*>(bc);
            if (!aleBc) continue;
            for (size_t faceId : aleBc->getFaceIds())
                faceBoundary[faceId] = aleBc;
        }
    }

    virtual void
    EvaluateDerivatives(const State<dim>* initialState,
                        State<dim>& deriv,
                        const double time,
                        const double dt) override {
        auto* rho = initialState->template getField<double>("density");
        auto* v   = initialState->template getField<Vector>("velocity");
        auto* u   = initialState->template getField<double>("specificInternalEnergy");

        auto* drhodt = deriv.template getField<double>("density");
        auto* dvdt   = deriv.template getField<Vector>("velocity");
        auto* dudt   = deriv.template getField<double>("specificInternalEnergy");

        auto* pressure   = this->nodeList->template getField<double>("pressure");
        auto* soundSpeed = this->nodeList->template getField<double>("soundSpeed");

        const auto& faces = mesh->getFaces();
        size_t numCells = mesh->getConnectivityMap().size();

        double local_dtmin = 1e30;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (size_t i = 0; i < numCells; ++i) {
            double rhoi = std::max(rho->getValue(i), 1e-12);
            Vector vi = v->getValue(i);
            double ui = std::max(u->getValue(i), 1e-12);
            double pi = pressure->getValue(i);
            double ci = soundSpeed->getValue(i);

            double net_rho_flux = 0.0;
            Vector net_mom_flux = Vector::zero();
            double net_E_flux = 0.0;
            double maxFaceArea = 0.0;

            for (size_t faceIdx : mesh->getFacesForCell(i)) {
                const Mesh::Face<dim>& f = faces[faceIdx];
                maxFaceArea = std::max(maxFaceArea, f.area);

                size_t other = (f.leftCell == i) ? f.rightCell : f.leftCell;

                double rhoOther, uOther;
                Vector vOther;
                if (other != Mesh::Face<dim>::boundaryCell) {
                    rhoOther = rho->getValue(other);
                    vOther   = v->getValue(other);
                    uOther   = u->getValue(other);
                } else {
                    ALEMeshBoundary<dim>* bc = faceBoundary[faceIdx];
                    if (bc) {
                        bc->ghostState(faceIdx, rhoi, vi, ui, rhoOther, vOther, uOther);
                    } else {
                        // No boundary registered for this face: default to a
                        // solid wall rather than silently leaking mass/energy.
                        double vn = vi.dot(f.normal);
                        vOther = vi - f.normal * (2.0 * vn);
                        rhoOther = rhoi;
                        uOther = ui;
                    }
                }

                double pOther, csOther;
                this->eos->setPressure(&pOther, &rhoOther, &uOther);
                this->eos->setSoundSpeed(&csOther, &rhoOther, &uOther);

                HLLFlux<dim> flux;
                if (f.leftCell == i) {
                    flux = this->computeFlux(rhoi, vi, ui, pi, ci,
                                              rhoOther, vOther, uOther, pOther, csOther,
                                              f.normal);
                    net_rho_flux -= flux.mass * f.area;
                    net_mom_flux -= flux.momentum * f.area;
                    net_E_flux   -= flux.energy * f.area;
                } else {
                    flux = this->computeFlux(rhoOther, vOther, uOther, pOther, csOther,
                                              rhoi, vi, ui, pi, ci,
                                              f.normal);
                    net_rho_flux += flux.mass * f.area;
                    net_mom_flux += flux.momentum * f.area;
                    net_E_flux   += flux.energy * f.area;
                }
            }

            double volume = mesh->cellVolume(i);
            double net_rho_rate = net_rho_flux / volume;
            Vector net_mom_rate = net_mom_flux / volume;
            double net_E_rate   = net_E_flux / volume;

            drhodt->setValue(i, net_rho_rate);
            Vector dvi = (net_mom_rate - vi * net_rho_rate) / rhoi;
            double dui = (net_E_rate - vi.dot(net_mom_rate) + 0.5 * vi.mag2() * net_rho_rate - ui * net_rho_rate) / rhoi;
            dvdt->setValue(i, dvi);
            dudt->setValue(i, dui);

            // Safety factor well below GridHydroBase's 0.1: a first-order (no
            // reconstruction) flux at a strong discontinuity, divided by a
            // small cell density, can produce a primitive-rate swing much
            // larger than the wave-speed-based CFL bound alone accounts for
            // (confirmed empirically on an 8x density-ratio Sod-style case --
            // 0.1 still diverged in a handful of steps, 0.02 did not).
            double L = volume / std::max(maxFaceArea, 1e-30);
            local_dtmin = std::min(local_dtmin, 0.02 * L / (ci + vi.magnitude() + 1e-30));
        }

        dtmin = local_dtmin;
    }

    virtual double
    EstimateTimestep() const override {
        return dtmin;
    }

    virtual void
    FinalizeStep(const State<dim>* finalState) override {
        auto* fdensity  = finalState->template getField<double>("density");
        auto* fvelocity = finalState->template getField<Vector>("velocity");
        auto* fu        = finalState->template getField<double>("specificInternalEnergy");

        auto* density  = this->nodeList->template getField<double>("density");
        auto* velocity = this->nodeList->template getField<Vector>("velocity");
        auto* u        = this->nodeList->template getField<double>("specificInternalEnergy");

        const double rho_floor = 1e-12;
        const double v_max = 1e3;
        const double u_max = 1e4;

        #pragma omp parallel for
        for (int i = 0; i < this->nodeList->size(); ++i) {
            double rhoi = std::max(fdensity->getValue(i), rho_floor);
            Vector vi = fvelocity->getValue(i);
            double ui = std::max(fu->getValue(i), rho_floor);

            if (rhoi <= 1.1 * rho_floor)
                vi = Vector::zero();
            else if (vi.magnitude() > v_max)
                vi = vi.normal() * v_max;

            if (ui > u_max)
                ui = u_max;

            density->setValue(i, rhoi);
            velocity->setValue(i, vi);
            u->setValue(i, ui);
        }

        EOSLookup();
    }

    virtual void
    EOSLookup() {
        NodeList* nodeList = this->nodeList;
        auto* rho = nodeList->getField<double>("density");
        auto* u = nodeList->getField<double>("specificInternalEnergy");
        auto* pressure = nodeList->getField<double>("pressure");
        auto* cs = nodeList->getField<double>("soundSpeed");
        this->eos->setPressure(pressure, rho, u);
        this->eos->setSoundSpeed(cs, rho, u);
    }

    // To be implemented by derived class -- given fully-reconstructed
    // primitive states on either side of a face (no reconstruction is done
    // here; see class comment), return the HLL-family flux across it.
    virtual HLLFlux<dim>
    computeFlux(double rhoL, const Vector& vL, double uL, double pL, double csL,
                double rhoR, const Vector& vR, double uR, double pR, double csR,
                const Vector& normal) const = 0;
};
