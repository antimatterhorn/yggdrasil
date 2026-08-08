// Copyright (C) 2026  Cody Raskin

#pragma once
#include <vector>
#include "../Mesh/grid.hh"
#include "../DataBase/nodeList.hh"
#include "../EOS/equationOfState.hh"
#include "../Physics/fluxObserver.hh"

namespace AMR {

// Accumulates one patch's face fluxes at the coarse-fine interface. Only faces
// registered via registerFace are kept; everything else the solver reports is
// ignored. Fluxes are stored as flux*area summed over derivative evaluations,
// with a per-slot count so the mean over an integrator's stages can be taken --
// for equal-weight stages (RK2) that mean times the coarse dt is exactly the
// material crossing the face over the step, including when a subcycled fine
// level contributes several substeps' worth of evaluations.
template <int dim>
class FluxRegister : public FluxObserver<dim> {
public:
    using Vector = Lin::Vector<dim>;

    FluxRegister(int numCells) : slotOf(numCells * dim * 2, -1) {}

    virtual ~FluxRegister() = default;

    int registerFace(int cell, int axis, bool plusSide) {
        const int key = keyOf(cell, axis, plusSide);
        if (slotOf[key] >= 0) return slotOf[key];
        const int slot = (int)cells.size();
        slotOf[key] = slot;
        cells.push_back(cell);
        plusSides.push_back(plusSide);
        massSum.push_back(0.0);
        momentumSum.push_back(Vector());
        energySum.push_back(0.0);
        counts.push_back(0);
        return slot;
    }

    void recordFlux(int cell, int axis, bool plusSide,
                    double massFlux, const Vector& momentumFlux, double energyFlux,
                    double area) override {
        const int key = keyOf(cell, axis, plusSide);
        if (key < 0 || key >= (int)slotOf.size()) return;
        const int slot = slotOf[key];
        if (slot < 0) return;
        massSum[slot]     += massFlux * area;
        momentumSum[slot] = momentumSum[slot] + momentumFlux * area;
        energySum[slot]   += energyFlux * area;
        counts[slot]      += 1;
    }

    void reset() {
        for (size_t s = 0; s < counts.size(); ++s) {
            massSum[s] = 0.0;
            momentumSum[s] = Vector();
            energySum[s] = 0.0;
            counts[s] = 0;
        }
    }

    int slotCount() const { return (int)cells.size(); }
    int cellOf(int slot) const { return cells[slot]; }
    bool plusSideOf(int slot) const { return plusSides[slot]; }
    double meanMass(int slot) const { return counts[slot] ? massSum[slot] / counts[slot] : 0.0; }
    Vector meanMomentum(int slot) const { return counts[slot] ? momentumSum[slot] * (1.0 / counts[slot]) : Vector(); }
    double meanEnergy(int slot) const { return counts[slot] ? energySum[slot] / counts[slot] : 0.0; }

private:
    int keyOf(int cell, int axis, bool plusSide) const { return (cell * dim + axis) * 2 + (plusSide ? 1 : 0); }

    std::vector<int> slotOf;
    std::vector<int> cells;
    std::vector<bool> plusSides;
    std::vector<double> massSum;
    std::vector<Vector> momentumSum;
    std::vector<double> energySum;
    std::vector<int> counts;
};

// Replaces the coarse solver's own interface flux with the sum of the fine
// fluxes across the same face, restoring conservation across a coarse-fine
// boundary. coarseSlots[k]/fineSlots[k] pair one coarse face with one of the
// fine faces covering it; the pairing is computed in AMRController.py.
template <int dim>
class Refluxer {
public:
    using Vector = Lin::Vector<dim>;

    Refluxer(NodeList* coarseNodes, Mesh::Grid<dim>* coarseGrid, EquationOfState* eos,
              FluxRegister<dim>* coarseRegister, FluxRegister<dim>* fineRegister,
              std::vector<int> coarseSlots, std::vector<int> fineSlots)
        : coarseNodes(coarseNodes), coarseGrid(coarseGrid), eos(eos),
          coarseRegister(coarseRegister), fineRegister(fineRegister),
          coarseSlots(std::move(coarseSlots)), fineSlots(std::move(fineSlots)) {}

    // Applies the correction for one coarse step, then clears both registers.
    void apply(double dt) {
        const int nSlots = coarseRegister->slotCount();
        std::vector<double> fineMass(nSlots, 0.0), fineEnergy(nSlots, 0.0);
        std::vector<Vector> fineMomentum(nSlots, Vector());
        std::vector<bool> touched(nSlots, false);

        for (size_t k = 0; k < coarseSlots.size(); ++k) {
            const int cs = coarseSlots[k];
            const int fs = fineSlots[k];
            fineMass[cs]     += fineRegister->meanMass(fs);
            fineMomentum[cs] = fineMomentum[cs] + fineRegister->meanMomentum(fs);
            fineEnergy[cs]   += fineRegister->meanEnergy(fs);
            touched[cs] = true;
        }

        auto* rho = coarseNodes->getFieldOrThrow<double>("density");
        auto* vel = coarseNodes->getFieldOrThrow<Vector>("velocity");
        auto* u   = coarseNodes->getFieldOrThrow<double>("specificInternalEnergy");

        for (int cs = 0; cs < nSlots; ++cs) {
            if (!touched[cs]) continue;
            // A +axis-side face enters the divergence with a minus sign, so its correction flips too.
            const double sign = coarseRegister->plusSideOf(cs) ? -1.0 : 1.0;
            const double dM   = sign * (fineMass[cs] - coarseRegister->meanMass(cs)) * dt;
            const Vector dP   = (fineMomentum[cs] - coarseRegister->meanMomentum(cs)) * (sign * dt);
            const double dE   = sign * (fineEnergy[cs] - coarseRegister->meanEnergy(cs)) * dt;

            const int c = coarseRegister->cellOf(cs);
            const double vol = coarseGrid->cellVolume(c);
            const Vector v0 = vel->getValue(c);
            const double m0 = rho->getValue(c) * vol;

            const double m = m0 + dM;
            if (m <= 0.0) continue;
            const Vector p = v0 * m0 + dP;
            const double e = m0 * (u->getValue(c) + 0.5 * v0.mag2()) + dE;

            const Vector v = p * (1.0 / m);
            rho->setValue(c, m / vol);
            vel->setValue(c, v);
            u->setValue(c, e / m - 0.5 * v.mag2());
        }

        auto* p  = coarseNodes->getFieldOrThrow<double>("pressure");
        auto* cs_ = coarseNodes->getFieldOrThrow<double>("soundSpeed");
        eos->setPressure(p, rho, u);
        eos->setSoundSpeed(cs_, rho, u);

        coarseRegister->reset();
        fineRegister->reset();
    }

private:
    NodeList* coarseNodes;
    Mesh::Grid<dim>* coarseGrid;
    EquationOfState* eos;
    FluxRegister<dim>* coarseRegister;
    FluxRegister<dim>* fineRegister;
    std::vector<int> coarseSlots;
    std::vector<int> fineSlots;
};

}
