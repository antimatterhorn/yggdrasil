// Copyright (C) 2026  Cody Raskin

#pragma once
#include "hydro.hh"
#include "../Mesh/grid.hh"
#include <iostream>
#include <unordered_set>

// Forward declaration for return type of computeFlux
template<int dim>
struct HLLFlux;

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
            if (!grid->onBoundary(i) && obstacleSet.count(i) == 0)
                insideIds.push_back(i);
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

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int h = 0; h < insideIds.size(); ++h) {
            int i = insideIds[h];
            Vector vi = v->getValue(i);
            double rhoi = rho->getValue(i);
            double ui = u->getValue(i);
            rhoi = std::max(rhoi, 1e-12);
            ui   = std::max(ui, 1e-12);
            double ci = soundSpeed->getValue(i);

            double ei = ui + 0.5 * vi.mag2();

            if (rhoi > 1e10 || ui > 1e10 || vi.magnitude() > 1e5) {
                std::cerr << "FATAL: Exploding state at cell " << i
                        << " rho=" << rhoi << " u=" << ui << " v=" << vi.toString() << std::endl;
                std::exit(EXIT_FAILURE);
            }

            Vector momi = vi * rhoi;
            double Ei = rhoi * ei;

            double net_rho_flux = 0.0;
            Vector net_mom_flux = Vector::zero();
            double net_E_flux = 0.0;

            auto neighbors = grid->getNeighboringCells(i);

            for (int k = 0; k < dim; ++k) {
                int jL = neighbors[2 * k];
                int jR = neighbors[2 * k + 1];

                // Check validity of jL and jR
                for (int q : {jL, jR, i}) {
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


                auto flux_L = this->computeFlux(jL, i, k, *rho, *v, *u, *pressure, *soundSpeed);
                auto flux_R = this->computeFlux(i, jR, k, *rho, *v, *u, *pressure, *soundSpeed);

                double dx = grid->spacing(k);
                net_rho_flux += (flux_L.mass - flux_R.mass) / dx;
                net_mom_flux += (flux_L.momentum - flux_R.momentum) / dx;
                net_E_flux   += (flux_L.energy - flux_R.energy) / dx;
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
        const double u_max = 1e4;

        #pragma omp parallel for
        for (int i = 0; i < this->nodeList->size(); ++i) {
            double rhoi = std::max(fdensity->getValue(i), rho_floor);
            Vector vi = fvelocity->getValue(i);
            double ui = std::max(fu->getValue(i), rho_floor);

            // Clamp velocity if rho is near floor
            if (rhoi <= 1.1 * rho_floor)
                vi = Vector::zero();
            else if (vi.magnitude() > v_max)
                vi = vi.normal() * v_max;

            // Clamp energy if absurd
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
