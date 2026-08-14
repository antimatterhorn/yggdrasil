// Copyright (C) 2026  Cody Raskin

#pragma once
#include "hydro.hh"
#include "../Mesh/grid.hh"
#include "../Boundaries/reflectingGridBoundary.cc"
#include <iostream>
#include <memory>
#include <unordered_set>

#include "HLL.cc"   // HLLFlux, the return type of computeFlux

template<int dim>
class GridHydroBase : public Hydro<dim> {
protected:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    Mesh::Grid<dim>* grid;
    std::vector<int> insideIds;
    double dxmin = 1e30;
    mutable double dtmin = 1e30;

    // A face is identified by the cell on its low side along the sweep axis, so cell i's
    // high face is face i and its low face is face neighbors[2k]. faceLowIds[k] lists
    // every face along axis k that some interior cell needs; faceFlux holds one axis'
    // worth of solved fluxes, indexed the same way and reused across axes.
    std::array<std::vector<int>, dim> faceLowIds;
    std::vector<HLLFlux<dim>> faceFlux;

    // r=0 axis (symmetry) boundary for RZ grids, owned here because the solver
    // installs it rather than the user. Null for Cartesian.
    std::unique_ptr<ReflectingGridBoundary<dim>> axisBoundary;

public:
    GridHydroBase(NodeList* nodeList,
                  PhysicalConstants& constants,
                  EquationOfState* eos,
                  Mesh::Grid<dim>* grid)
        : Hydro<dim>(nodeList, constants, eos), grid(grid) {
        
        grid->assignPositions(nodeList);
        State<dim>* state = &this->state;

        auto* v   = nodeList->getField<Vector>("velocity");
        auto* rho = nodeList->getField<double>("density");
        auto* u   = nodeList->getField<double>("specificInternalEnergy");

        state->template addField<Vector>(v);
        state->template addField<double>(rho);
        state->template addField<double>(u);

        for (int i = 0; i < dim; ++i)
            dxmin = std::min(dxmin, grid->spacing(i));
    }

    virtual ~GridHydroBase() {}

    virtual void
    ZeroTimeInitialize() override {
        EOSLookup();
        State<dim> state = this->state;
        NodeList* nodeList = this->nodeList;
        this->UpdateState();

        // The r=0 axis is a symmetry line, i.e. a reflecting condition on the
        // r-min ("left", face 0) boundary. Reject a user boundary already on that
        // face, then append ours so it applies last.
        if (grid->geometry() == Mesh::Geometry::CylindricalRZ && !axisBoundary) {
            for (auto* bc : this->boundaries) {
                auto* gb = dynamic_cast<GridBoundary<dim>*>(bc);
                if (gb && gb->isFaceActive(0))
                    throw std::runtime_error(
                        "GridHydroBase: in CylindricalRZ geometry the r=0 axis "
                        "(the 'left' face) is managed automatically; do not "
                        "assign a boundary to it.");
            }
            axisBoundary = std::make_unique<ReflectingGridBoundary<dim>>(grid);
            axisBoundary->setFaces({"left"});
            this->addBoundary(axisBoundary.get());
        }

        // InitializeBoundaries must run first so that any ReflectingGridBoundary
        // obstacles have their IDs deduped and their neighbor lists built before
        // we query getObstacleIds() below.
        this->InitializeBoundaries();

        // Collect interior obstacle cells from all registered boundaries and
        // exclude them from the flux update loop, mirroring WaveEquation's
        // obstacle handling — otherwise a "solid" cell would still have fluxes
        // computed into/out of it as if it were ordinary fluid. Excluded cells
        // are left for the owning boundary's own Apply (e.g. ReflectingGridBoundary's
        // interior Neumann fill) to keep finite each step.
        std::unordered_set<int> obstacleSet;
        for (auto* bc : this->boundaries)
            for (int id : bc->getObstacleIds())
                obstacleSet.insert(id);

        insideIds.clear();
        for (int i = 0; i < grid->size(); ++i)
            if (!grid->isGhost(i) && obstacleSet.count(i) == 0)
                insideIds.push_back(i);

        // Every face an interior cell touches, per axis. An interior cell is never on the
        // rim, so its high neighbour always exists; a low neighbour pulled in here is at
        // worst a rim cell, whose high neighbour is the interior cell that named it. So
        // every face listed has a valid pair.
        faceFlux.assign(grid->size(), HLLFlux<dim>{});
        std::vector<char> needed(grid->size());
        for (int k = 0; k < dim; ++k) {
            std::fill(needed.begin(), needed.end(), 0);
            for (int i : insideIds) {
                needed[i] = 1;
                const int jL = grid->getNeighboringCells(i)[2 * k];
                if (jL >= 0) needed[jL] = 1;
            }
            faceLowIds[k].clear();
            for (int j = 0; j < grid->size(); ++j)
                if (needed[j]) faceLowIds[k].push_back(j);
        }
    }

    virtual void 
    EvaluateDerivatives(const State<dim>* initialState,
                                     State<dim>& deriv,
                                     const double time,
                                     const double dt) override {
        NodeList* nodeList = this->nodeList;

        auto* rho = initialState->template getField<double>("density");
        auto* v   = initialState->template getField<Vector>("velocity");
        auto* u   = initialState->template getField<double>("specificInternalEnergy");

        auto* drhodt = deriv.template getField<double>("density");
        auto* dvdt   = deriv.template getField<Vector>("velocity");
        auto* dudt   = deriv.template getField<double>("specificInternalEnergy");

        auto* pressure   = nodeList->getField<double>("pressure");
        auto* soundSpeed = nodeList->getField<double>("soundSpeed");

        double local_dtmin = 1e30;
        const int nInside = (int)insideIds.size();

        // State sanity, once per interior cell and its neighbours rather than once per
        // axis. A rim cell is covered as some interior cell's neighbour.
        #pragma omp parallel for
        for (int h = 0; h < nInside; ++h) {
            const int i = insideIds[h];
            if (rho->getValue(i) > 1e10 || u->getValue(i) > 1e10 ||
                v->getValue(i).mag2() > 1e10) {
                std::cerr << "FATAL: Exploding state at cell " << i
                        << " rho=" << rho->getValue(i) << " u=" << u->getValue(i)
                        << " v=" << v->getValue(i).toString() << std::endl;
                std::exit(EXIT_FAILURE);
            }
            const auto neighbors = grid->getNeighboringCells(i);
            for (int k = 0; k < 2 * dim; ++k) {
                const int q = neighbors[k];
                if (q < 0) continue;
                if (!std::isfinite(rho->getValue(q)) ||
                    !std::isfinite(u->getValue(q)) ||
                    !std::isfinite(v->getValue(q)[0])) {
                    std::cerr << "BAD INPUT at cell " << q << ": "
                            << "rho=" << rho->getValue(q)
                            << ", u=" << u->getValue(q)
                            << ", vx=" << v->getValue(q)[0] << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
        }

        // Net flux accumulators, kept in the derivative fields until the conversion
        // pass below.
        #pragma omp parallel for
        for (int h = 0; h < nInside; ++h) {
            const int i = insideIds[h];
            drhodt->setValue(i, 0.0);
            dvdt->setValue(i, Vector::zero());
            dudt->setValue(i, 0.0);
        }

        for (int k = 0; k < dim; ++k) {
            // Solve each face once. The old form evaluated cell i's high face and cell
            // jR's low face separately, which are the same face with the same arguments.
            const auto& faces = faceLowIds[k];
            #pragma omp parallel for
            for (int h = 0; h < (int)faces.size(); ++h) {
                const int j  = faces[h];
                const int jR = grid->getNeighboringCells(j)[2 * k + 1];
                faceFlux[j] = this->computeFlux(j, jR, k, *rho, *v, *u, *pressure, *soundSpeed);
            }

            // Area-weighted finite-volume divergence: sum(flux * face area) / cell
            // volume. With the RZ metric this is the cylindrical divergence
            // (1/r) d(r F)/dr + dF/dz.
            #pragma omp parallel for
            for (int h = 0; h < nInside; ++h) {
                const int i  = insideIds[h];
                const auto neighbors = grid->getNeighboringCells(i);
                const int jL = neighbors[2 * k];
                const int jR = neighbors[2 * k + 1];

                const double Vi  = grid->cellVolume(i);
                const double A_L = grid->faceArea(jL, i, k);
                const double A_R = grid->faceArea(i, jR, k);

                const HLLFlux<dim>& flux_L = faceFlux[jL];
                const HLLFlux<dim>& flux_R = faceFlux[i];

                drhodt->setValue(i, drhodt->getValue(i) +
                                 (flux_L.mass * A_L - flux_R.mass * A_R) / Vi);
                dvdt->setValue(i, dvdt->getValue(i) +
                               (flux_L.momentum * A_L - flux_R.momentum * A_R) / Vi);
                dudt->setValue(i, dudt->getValue(i) +
                               (flux_L.energy * A_L - flux_R.energy * A_R) / Vi);
            }
        }

        // Convert the accumulated net fluxes into the conserved-variable derivatives.
        #pragma omp parallel for reduction(min:local_dtmin)
        for (int h = 0; h < nInside; ++h) {
            const int i = insideIds[h];
            const Vector vi = v->getValue(i);
            const double rhoi = std::max(rho->getValue(i), 1e-12);
            const double ui   = std::max(u->getValue(i), 1e-12);
            const double ci   = soundSpeed->getValue(i);

            double net_rho_flux = drhodt->getValue(i);
            Vector net_mom_flux = dvdt->getValue(i);
            double net_E_flux   = dudt->getValue(i);

            // Cylindrical geometric source: the radius-weighted divergence above
            // puts an extra p/r on radial momentum ((1/r) d(r*p)/dr = dp/dr +
            // p/r), so adding p/r back leaves only the physical dp/dr. Written
            // p*(A_Rr - A_Lr)/Vi (= p/r) to reuse the same radial face areas and
            // volume as the flux, so a uniform gas stays exactly at rest.
            if (grid->geometry() == Mesh::Geometry::CylindricalRZ) {
                const auto neighbors = grid->getNeighboringCells(i);
                const double A_Lr = grid->faceArea(neighbors[0], i, 0);
                const double A_Rr = grid->faceArea(i, neighbors[1], 0);
                net_mom_flux[0] += pressure->getValue(i) * (A_Rr - A_Lr) / grid->cellVolume(i);
            }

            drhodt->setValue(i, net_rho_flux);
            Vector dvi = (net_mom_flux - vi * net_rho_flux) / rhoi;
            double dui = (net_E_flux - vi.dot(net_mom_flux) + 0.5 * vi.mag2() * net_rho_flux - ui * net_rho_flux) / rhoi;

            dvdt->setValue(i, dvi);
            dudt->setValue(i, dui);

            local_dtmin = std::min(local_dtmin, 0.1 * dxmin / (ci+vi.magnitude()));
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
        const double u_max = 1e6;

        #pragma omp parallel for
        for (int i = 0; i < this->nodeList->size(); ++i) {
            double rhoi = std::max(fdensity->getValue(i), rho_floor);
            Vector vi = fvelocity->getValue(i);
            double ui = std::max(fu->getValue(i), rho_floor);

            // Clamp velocity if rho is near floor
            if (rhoi <= 1.1 * rho_floor)
                vi = Vector::zero();
            else if (vi.mag2() > v_max * v_max)
                vi = vi.normal() * v_max;

            // Clamp energy if absurd
            if (ui > u_max) {
                std::cout << "GridHydroBase: energy cap triggered at cell " << i
                          << " (u=" << ui << " > u_max=" << u_max << ")" << std::endl;
                ui = u_max;
            }

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

    virtual double 
    getCell(int i, int j, const std::string& fieldName = "pressure") const {
        int idx = grid->index(i, j, 0);
        ScalarField* field = this->nodeList->template getField<double>(fieldName);
        return field->getValue(idx);
    }

    virtual double 
    getCellComponent(int i, int j, int component, const std::string& fieldName) const {
        int idx = grid->index(i, j, 0);
        auto* field = this->nodeList->template getField<Vector>(fieldName);
        return field->getValue(idx)[component];
    }

    // To be implemented by derived class
    virtual HLLFlux<dim> 
    computeFlux(int iL, int iR, int axis,
                const Field<double>& rho,
                const Field<Vector>& v,
                const Field<double>& u,
                const Field<double>& p,
                const Field<double>& cs) const = 0;
};
