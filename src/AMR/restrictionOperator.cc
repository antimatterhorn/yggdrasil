// Copyright (C) 2026  Cody Raskin

#pragma once
#include <vector>
#include <algorithm>
#include "../Mesh/grid.hh"
#include "../DataBase/nodeList.hh"
#include "../EOS/equationOfState.hh"

namespace AMR {

// Averages a fine patch's interior cells down onto the coarse cells they cover.
// Sums extensive mass/momentum/total energy and divides back out rather than
// averaging density/velocity/sie directly, which would not conserve momentum
// or energy. fineIds[k] contributes to coarseIds[k]; the pairing (only fully
// covered coarse cells) is computed in AMRController.py.
template <int dim>
class RestrictionOperator {
public:
    using Vector = Lin::Vector<dim>;

    RestrictionOperator(NodeList* fineNodes, NodeList* coarseNodes,
                         Mesh::Grid<dim>* fineGrid, Mesh::Grid<dim>* coarseGrid,
                         EquationOfState* eos,
                         std::vector<int> fineIds, std::vector<int> coarseIds)
        : fineNodes(fineNodes), coarseNodes(coarseNodes),
          fineGrid(fineGrid), coarseGrid(coarseGrid), eos(eos),
          fineIds(std::move(fineIds)), coarseIds(std::move(coarseIds)) {
        distinctCoarseIds = this->coarseIds;
        std::sort(distinctCoarseIds.begin(), distinctCoarseIds.end());
        distinctCoarseIds.erase(std::unique(distinctCoarseIds.begin(), distinctCoarseIds.end()),
                                 distinctCoarseIds.end());
    }

    void apply() {
        auto* fRho = fineNodes->template getFieldOrThrow<double>("density");
        auto* fVel = fineNodes->template getFieldOrThrow<Vector>("velocity");
        auto* fU   = fineNodes->template getFieldOrThrow<double>("specificInternalEnergy");

        auto* cRho = coarseNodes->template getFieldOrThrow<double>("density");
        auto* cVel = coarseNodes->template getFieldOrThrow<Vector>("velocity");
        auto* cU   = coarseNodes->template getFieldOrThrow<double>("specificInternalEnergy");

        const int nCoarse = (int)coarseNodes->size();
        std::vector<double> mass(nCoarse, 0.0), energy(nCoarse, 0.0);
        std::vector<Vector> momentum(nCoarse, Vector());

        for (size_t k = 0; k < fineIds.size(); ++k) {
            const int f = fineIds[k];
            const int c = coarseIds[k];
            const double vol = fineGrid->cellVolume(f);
            const double m   = fRho->getValue(f) * vol;
            const Vector v   = fVel->getValue(f);
            mass[c]     += m;
            momentum[c] = momentum[c] + v * m;
            energy[c]   += m * (fU->getValue(f) + 0.5 * v.mag2());
        }

        for (const int c : distinctCoarseIds) {
            if (mass[c] <= 0.0) continue;
            const Vector v = momentum[c] * (1.0 / mass[c]);
            cRho->setValue(c, mass[c] / coarseGrid->cellVolume(c));
            cVel->setValue(c, v);
            cU->setValue(c, energy[c] / mass[c] - 0.5 * v.mag2());
        }

        // Pressure/soundSpeed would otherwise stay stale for every restricted cell.
        auto* cP  = coarseNodes->template getFieldOrThrow<double>("pressure");
        auto* cCs = coarseNodes->template getFieldOrThrow<double>("soundSpeed");
        eos->setPressure(cP, cRho, cU);
        eos->setSoundSpeed(cCs, cRho, cU);
    }

private:
    NodeList* fineNodes;
    NodeList* coarseNodes;
    Mesh::Grid<dim>* fineGrid;
    Mesh::Grid<dim>* coarseGrid;
    EquationOfState* eos;
    std::vector<int> fineIds;
    std::vector<int> coarseIds;
    std::vector<int> distinctCoarseIds;
};

}
