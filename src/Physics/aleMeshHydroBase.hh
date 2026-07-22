// Copyright (C) 2026  Cody Raskin

#pragma once
#include "hydro.hh"
#include "../Mesh/alemesh.hh"
#include "../Boundaries/aleMeshBoundary.hh"
#include "../Boundaries/reflectingALEMeshBoundary.cc"
#include <iostream>
#include <algorithm>
#include <vector>
#include <memory>
#include <stdexcept>

// Forward declaration for return type of computeFlux (definition in HLL.cc,
// included by concrete solvers -- same pattern as GridHydroBase).
template<int dim>
struct HLLFlux;

// Finite-volume ALE hydro on Mesh::ALEMesh, walking Face objects instead of
// Grid's axis-aligned neighbor pairs. One NodeList entry per ALEMesh cell
// (cell-centered, positions taken from mesh->cellCentroid), same convention
// as GridHydroBase.
//
// Tracks EXTENSIVE per-cell quantities (mass, momentum, total energy) as the
// integrated State fields, not intensive density/velocity/specific-internal-
// energy. This is what makes real mesh motion possible within the existing
// RK-integrator/State<dim> machinery: State's arithmetic (+=, *=scalar) is a
// pure per-field linear combination with no notion of geometry, so it can
// only correctly integrate a quantity whose ODE has no explicit volume term.
// Extensive quantities have exactly that property -- d(mass)/dt is a pure
// flux divergence regardless of whether the cell's volume is changing,
// *provided* the flux used is the ALE (moving-face) flux (see HLL.cc's
// computeHLLCFluxFromStatesALE / computeHLLEFluxFromStatesALE). Density/
// velocity/specificInternalEnergy are derived quantities (extensive /
// mesh->cellVolume(i)), recomputed each EvaluateDerivatives call and again
// in FinalizeStep.
//
// Mesh velocity, if any, is user-prescribed via an optional node-centered
// NodeList (one entry per mesh node, not per cell) holding a "velocity"
// Vector field -- zero everywhere for pure Eulerian (the default, matching
// prior behavior exactly, since S=0 makes the ALE flux reduce algebraically
// to the static flux), copied from material velocity for Lagrangian, or
// anything in between. No automatic rezoning; see root CLAUDE.md. The mesh
// itself only moves once per full timestep (in FinalizeStep, after the RK
// sub-stepping above has completed) -- geometry is fixed for the duration of
// one Controller.Step()'s integration, same as it always effectively was.
template<int dim>
class ALEMeshHydroBase : public Hydro<dim> {
protected:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    Mesh::ALEMesh<dim>* mesh;
    mutable double dtmin = 1e30;

    NodeList* nodeVelocities;
    VectorField* meshVelocity = nullptr;

    std::vector<ALEMeshBoundary<dim>*> faceBoundary;
    std::unique_ptr<ReflectingALEMeshBoundary<dim>> axisBoundary;

public:
    ALEMeshHydroBase(NodeList* nodeList,
                      PhysicalConstants& constants,
                      EquationOfState* eos,
                      Mesh::ALEMesh<dim>* mesh,
                      NodeList* nodeVelocities = nullptr)
        : Hydro<dim>(nodeList, constants, eos), mesh(mesh), nodeVelocities(nodeVelocities) {

        size_t numCells = mesh->getConnectivityMap().size();
        auto* pos = nodeList->getField<Vector>("position");
        for (size_t i = 0; i < numCells; ++i)
            pos->setValue(i, mesh->cellCentroid(i));

        if (nodeVelocities) {
            if (nodeVelocities->getField<Vector>("velocity") == nullptr)
                nodeVelocities->insertField<Vector>("velocity");
            meshVelocity = nodeVelocities->getField<Vector>("velocity");
        }

        // Extensive per-cell conserved quantities -- the actual INTEGRATE
        // fields. density/velocity/specificInternalEnergy (already enrolled
        // by Hydro<dim>/VerifyFields) become derived-only: still present on
        // the NodeList for EOS lookups, boundary ghostState calls, and IO,
        // but no longer tracked in `state`.
        this->template EnrollFields<double>({"cellMass", "cellEnergy"});
        this->template EnrollFields<Vector>({"cellMomentum"});
        this->template EnrollStateFields<double>({"cellMass", "cellEnergy"});
        this->template EnrollStateFields<Vector>({"cellMomentum"});
    }

    virtual ~ALEMeshHydroBase() {}

    virtual void
    ZeroTimeInitialize() override {
        EOSLookup();

        if (mesh->getFaces().empty())
            mesh->computeFaces();

        if (mesh->geometry() == Mesh::Geometry::CylindricalRZ) {
            const double tiny = 1e-12;
            const auto& nodePositions = mesh->getNodes();

            for (const auto& p : nodePositions)
                if (p.x() < -tiny)
                    throw std::runtime_error(
                        "ALEMeshHydroBase: CylindricalRZ geometry requires every "
                        "mesh node to have radius (x()) >= 0.");

            // The r=0 axis is a symmetry line: any boundary face whose two
            // endpoint nodes both lie at r=0. Reject a user boundary already
            // covering one of those faces, then install ours so it applies
            // last -- mirrors GridHydroBase's "left" face guard exactly.
            if (!axisBoundary) {
                std::vector<size_t> axisFaces;
                const auto& faces = mesh->getFaces();
                for (size_t faceId = 0; faceId < faces.size(); ++faceId) {
                    const auto& f = faces[faceId];
                    if (!f.isBoundary()) continue;
                    bool onAxis = true;
                    for (size_t nid : f.nodeIndices)
                        if (nodePositions[nid].x() > tiny) onAxis = false;
                    if (onAxis) axisFaces.push_back(faceId);
                }

                for (auto* bc : this->boundaries) {
                    auto* aleBc = dynamic_cast<ALEMeshBoundary<dim>*>(bc);
                    if (!aleBc) continue;
                    for (size_t faceId : aleBc->getFaceIds())
                        if (std::find(axisFaces.begin(), axisFaces.end(), faceId) != axisFaces.end())
                            throw std::runtime_error(
                                "ALEMeshHydroBase: in CylindricalRZ geometry the r=0 "
                                "axis faces are managed automatically; do not assign "
                                "a boundary to them.");
                }

                if (!axisFaces.empty()) {
                    axisBoundary = std::make_unique<ReflectingALEMeshBoundary<dim>>(mesh);
                    axisBoundary->setFaces(axisFaces);
                    this->addBoundary(axisBoundary.get());
                }
            }
        }

        // Seed the extensive fields from whatever density/velocity/specific-
        // internal-energy the caller set as the initial condition.
        auto* rho  = this->nodeList->template getField<double>("density");
        auto* v    = this->nodeList->template getField<Vector>("velocity");
        auto* u    = this->nodeList->template getField<double>("specificInternalEnergy");
        auto* mass = this->nodeList->template getField<double>("cellMass");
        auto* mom  = this->nodeList->template getField<Vector>("cellMomentum");
        auto* E    = this->nodeList->template getField<double>("cellEnergy");

        size_t numCells = mesh->getConnectivityMap().size();
        for (size_t i = 0; i < numCells; ++i) {
            double rhoi = rho->getValue(i);
            Vector vi = v->getValue(i);
            double ui = u->getValue(i);
            double Vi = mesh->cellVolume(i);
            mass->setValue(i, rhoi * Vi);
            mom->setValue(i, (rhoi * Vi) * vi);
            E->setValue(i, rhoi * Vi * (ui + 0.5 * vi.mag2()));
        }

        this->UpdateState();
        this->InitializeBoundaries();

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
        // Recorded so FinalizeStep can move the mesh by exactly this step's
        // dt: every concrete integrator's *last* EvaluateDerivatives call
        // within one Integrate() passes the full step dt (RK4's k4, RK2's
        // k2), and FinalizeStep always runs immediately after Integrate()
        // completes, so by the time FinalizeStep reads this->lastDt it holds
        // this step's exact dt, not an approximation.
        this->lastDt = dt;

        auto* mass = initialState->template getField<double>("cellMass");
        auto* mom  = initialState->template getField<Vector>("cellMomentum");
        auto* E    = initialState->template getField<double>("cellEnergy");

        auto* dmassdt = deriv.template getField<double>("cellMass");
        auto* dmomdt  = deriv.template getField<Vector>("cellMomentum");
        auto* dEdt    = deriv.template getField<double>("cellEnergy");

        // "Self" pressure/soundSpeed come from the NodeList fields (last set
        // by EOSLookup at the previous FinalizeStep), not recomputed per RK
        // sub-stage like the freshly-derived rho/v/u below.
        auto* pressure   = this->nodeList->template getField<double>("pressure");
        auto* soundSpeed = this->nodeList->template getField<double>("soundSpeed");

        const auto& faces = mesh->getFaces();
        size_t numCells = mesh->getConnectivityMap().size();

        // Derive this sub-stage's intensive primitive state (rho, v, u) from
        // the extensive quantities once per cell -- every face touching a
        // cell needs it, so compute it up front rather than re-deriving it
        // (an extensive/volume division) on every face visit.
        std::vector<double> rhoArr(numCells), uArr(numCells);
        std::vector<Vector> vArr(numCells);

        #pragma omp parallel for
        for (size_t i = 0; i < numCells; ++i) {
            double Vi = mesh->cellVolume(i);
            double mi = std::max(mass->getValue(i), 1e-12);
            Vector momi = mom->getValue(i);
            double Ei = E->getValue(i);

            double rhoi = mi / Vi;
            Vector vi = momi / mi;
            double ui = std::max(Ei / mi - 0.5 * vi.mag2(), 1e-12);

            rhoArr[i] = rhoi;
            vArr[i] = vi;
            uArr[i] = ui;
        }

        double local_dtmin = 1e30;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (size_t i = 0; i < numCells; ++i) {
            double rhoi = std::max(rhoArr[i], 1e-12);
            Vector vi = vArr[i];
            double ui = std::max(uArr[i], 1e-12);
            double pi = pressure->getValue(i);
            double ci = soundSpeed->getValue(i);

            double net_mass_flux = 0.0;
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
                    rhoOther = rhoArr[other];
                    vOther   = vArr[other];
                    uOther   = uArr[other];
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

                // Face normal speed S = w.n (the mesh's own prescribed
                // motion), averaged over the face's endpoint nodes -- zero
                // everywhere unless a node-velocity NodeList was given.
                double S = 0.0;
                if (meshVelocity) {
                    Vector wSum = Vector::zero();
                    for (size_t nid : f.nodeIndices)
                        wSum = wSum + meshVelocity->getValue(nid);
                    S = (wSum / static_cast<double>(f.nodeIndices.size())).dot(f.normal);
                }

                HLLFlux<dim> flux;
                if (f.leftCell == i) {
                    flux = this->computeFlux(rhoi, vi, ui, pi, ci,
                                              rhoOther, vOther, uOther, pOther, csOther,
                                              f.normal, S);
                    net_mass_flux -= flux.mass * f.area;
                    net_mom_flux  -= flux.momentum * f.area;
                    net_E_flux    -= flux.energy * f.area;
                } else {
                    flux = this->computeFlux(rhoOther, vOther, uOther, pOther, csOther,
                                              rhoi, vi, ui, pi, ci,
                                              f.normal, S);
                    net_mass_flux += flux.mass * f.area;
                    net_mom_flux  += flux.momentum * f.area;
                    net_E_flux    += flux.energy * f.area;
                }
            }

            // Cylindrical geometric source: the radius-weighted face areas
            // above make the divergence pick up an extra p/r on radial
            // momentum ((1/r) d(r p)/dr = dp/dr + p/r), so adding
            // p_cell * A_cell,2D back (the extensive form of GridHydroBase's
            // intensive p/r term) leaves only the physical dp/dr. A_cell,2D
            // is recovered from the RZ-weighted cellVolume rather than
            // adding a separate mesh accessor for the plain area.
            if (mesh->geometry() == Mesh::Geometry::CylindricalRZ) {
                double ri = mesh->cellCentroid(i).x();
                net_mom_flux[0] += pi * mesh->cellVolume(i) / ri;
            }

            dmassdt->setValue(i, net_mass_flux);
            dmomdt->setValue(i, net_mom_flux);
            dEdt->setValue(i, net_E_flux);

            // Safety factor well below GridHydroBase's 0.1: a first-order (no
            // reconstruction) flux at a strong discontinuity, divided by a
            // small cell density, can produce a primitive-rate swing much
            // larger than the wave-speed-based CFL bound alone accounts for
            // (confirmed empirically on an 8x density-ratio Sod-style case --
            // 0.1 still diverged in a handful of steps, 0.02 did not).
            double volume = mesh->cellVolume(i);
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
        auto* fmass = finalState->template getField<double>("cellMass");
        auto* fmom  = finalState->template getField<Vector>("cellMomentum");
        auto* fE    = finalState->template getField<double>("cellEnergy");

        auto* mass = this->nodeList->template getField<double>("cellMass");
        auto* mom  = this->nodeList->template getField<Vector>("cellMomentum");
        auto* E    = this->nodeList->template getField<double>("cellEnergy");

        const double mass_floor = 1e-12;

        #pragma omp parallel for
        for (int i = 0; i < this->nodeList->size(); ++i) {
            mass->setValue(i, std::max(fmass->getValue(i), mass_floor));
            mom->setValue(i, fmom->getValue(i));
            E->setValue(i, fE->getValue(i));
        }

        // Move the mesh once per full step, not once per RK sub-stage --
        // geometry (face areas/normals, cell volumes) used throughout the RK
        // sub-stepping above was the geometry at the START of this step;
        // only now does it advance to match the end of the step.
        if (meshVelocity) {
            const auto& nodePositions = mesh->getNodes();
            for (size_t i = 0; i < nodePositions.size(); ++i)
                mesh->setNodePosition(i, nodePositions[i] + meshVelocity->getValue(i) * this->lastDt);
            mesh->updateGeometry();

            auto* pos = this->nodeList->template getField<Vector>("position");
            for (int i = 0; i < this->nodeList->size(); ++i)
                pos->setValue(i, mesh->cellCentroid(i));
        }

        // Derive density/velocity/specificInternalEnergy from the just-
        // integrated extensive fields using the mesh's *current* (just
        // updated, if it moved) cell volumes -- this division is the entire
        // "remap": no separate geometric projection needed, since the
        // extensive quantities were already correctly advanced through the
        // ALE flux above.
        auto* density  = this->nodeList->template getField<double>("density");
        auto* velocity = this->nodeList->template getField<Vector>("velocity");
        auto* u        = this->nodeList->template getField<double>("specificInternalEnergy");

        const double rho_floor = 1e-12;
        const double v_max = 1e3;
        const double u_max = 1e4;

        #pragma omp parallel for
        for (int i = 0; i < this->nodeList->size(); ++i) {
            double Vi = mesh->cellVolume(i);
            double mi = std::max(mass->getValue(i), mass_floor);
            double rhoi = std::max(mi / Vi, rho_floor);
            Vector vi = mom->getValue(i) / mi;
            double ui = std::max(E->getValue(i) / mi - 0.5 * vi.mag2(), rho_floor);

            if (rhoi <= 1.1 * rho_floor)
                vi = Vector::zero();
            else if (vi.magnitude() > v_max)
                vi = vi.normal() * v_max;

            if (ui > u_max)
                ui = u_max;

            density->setValue(i, rhoi);
            velocity->setValue(i, vi);
            u->setValue(i, ui);

            // Keep the extensive fields consistent with any flooring/
            // clamping just applied, so the next step's derivative
            // evaluation starts from the same (possibly clamped) state.
            mass->setValue(i, rhoi * Vi);
            mom->setValue(i, (rhoi * Vi) * vi);
            E->setValue(i, rhoi * Vi * (ui + 0.5 * vi.mag2()));
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
    // here; see class comment) and the face's own normal speed S = w.n
    // (zero for a static mesh), return the HLL-family ALE flux across it.
    virtual HLLFlux<dim>
    computeFlux(double rhoL, const Vector& vL, double uL, double pL, double csL,
                double rhoR, const Vector& vR, double uR, double pR, double csR,
                const Vector& normal, double faceVelocity) const = 0;
};
