// Copyright (C) 2026  Cody Raskin

// GridHydroHLLC.hh
#pragma once
#include "gridHydroBase.hh"
#include "HLL.cc"

template<int dim>
class GridHydroHLLC : public GridHydroBase<dim> {
public:
    using typename GridHydroBase<dim>::Vector;
    using typename GridHydroBase<dim>::VectorField;
    using typename GridHydroBase<dim>::ScalarField;

    GridHydroHLLC(NodeList* nodeList, PhysicalConstants& constants,
                  EquationOfState* eos, Mesh::Grid<dim>* grid)
        : GridHydroBase<dim>(nodeList, constants, eos, grid) {}

    virtual std::string name() const override { return "GridHydroHLLC"; }

    virtual
    HLLFlux<dim> 
    computeFlux(int iL, int iR, int axis,
        const Field<double>& rho, const Field<Vector>& v,
        const Field<double>& u, const Field<double>& p,
        const Field<double>& cs) const override {

        using Vector = Lin::Vector<dim>;
        using Lin::minmod;
        using Lin::minmodVec;

        auto neighborsL = this->grid->getNeighboringCells(iL);
        auto neighborsR = this->grid->getNeighboringCells(iR);

        // Near a domain edge the second-layer neighbor may not exist (-1);
        // fall back to the immediate neighbor so the corresponding minmod
        // difference collapses to zero instead of reading an invalid index.
        int iLL = neighborsL[2 * axis];
        int iRR = neighborsR[2 * axis + 1];
        const bool atEdge = (iLL < 0 || iRR < 0);
        if (iLL < 0) iLL = iL;
        if (iRR < 0) iRR = iR;

        // Center values
        double rhoL0 = rho.getValue(iL);
        double rhoR0 = rho.getValue(iR);
        double uL0 = u.getValue(iL);
        double uR0 = u.getValue(iR);
        Vector vL0 = v.getValue(iL);
        Vector vR0 = v.getValue(iR);

        // Slopes (minmod-limited)
        double srhoL = minmod(rho.getValue(iL) - rho.getValue(iLL), rho.getValue(iR) - rho.getValue(iL));
        double srhoR = minmod(rho.getValue(iR) - rho.getValue(iL), rho.getValue(iRR) - rho.getValue(iR));

        double suL = minmod(u.getValue(iL) - u.getValue(iLL), u.getValue(iR) - u.getValue(iL));
        double suR = minmod(u.getValue(iR) - u.getValue(iL), u.getValue(iRR) - u.getValue(iR));

        Vector svL = minmodVec(v.getValue(iL) - v.getValue(iLL), v.getValue(iR) - v.getValue(iL));
        Vector svR = minmodVec(v.getValue(iR) - v.getValue(iL), v.getValue(iRR) - v.getValue(iR));

        // Optional: clamp slope values to avoid overshoots
        const double slopeLimiterMax = 10.0;
        srhoL = std::clamp(srhoL, -slopeLimiterMax, slopeLimiterMax);
        srhoR = std::clamp(srhoR, -slopeLimiterMax, slopeLimiterMax);
        suL = std::clamp(suL, -slopeLimiterMax, slopeLimiterMax);
        suR = std::clamp(suR, -slopeLimiterMax, slopeLimiterMax);
        for (int d = 0; d < dim; ++d) {
            svL[d] = std::clamp(svL[d], -slopeLimiterMax, slopeLimiterMax);
            svR[d] = std::clamp(svR[d], -slopeLimiterMax, slopeLimiterMax);
        }

        // A one-sided slope at a domain edge breaks the mirror symmetry a wall's flux cancellation needs.
        if (atEdge) {
            srhoL = srhoR = 0.0;
            suL   = suR   = 0.0;
            svL   = svR   = Vector::zero();
        }

        // Reconstruct states at interface
        double rhoL = rhoL0 + 0.5 * srhoL;
        double rhoR = rhoR0 - 0.5 * srhoR;
        double uL   = uL0   + 0.5 * suL;
        double uR   = uR0   - 0.5 * suR;
        Vector vL   = vL0   + 0.5 * svL;
        Vector vR   = vR0   - 0.5 * svR;

        // Clamp to physical floors
        rhoL = std::max(rhoL, 1e-12);
        rhoR = std::max(rhoR, 1e-12);
        uL = std::clamp(uL, 1e-12, 1e4);
        uR = std::clamp(uR, 1e-12, 1e4);

        // Clamp reconstructed velocities. Tested squared so the sqrt is paid only by the
        // rare vector that actually needs clamping.
        const double vmax = 1e3;
        if (vL.mag2() > vmax * vmax) vL = vL.unit() * vmax;
        if (vR.mag2() > vmax * vmax) vR = vR.unit() * vmax;

        // Pressure and sound speed from the reconstructed (rho, u) via the EOS,
        // so the interface state is thermodynamically consistent (the passed-in
        // cell-centered p/cs do not match the reconstructed rho/u). p/cs unused.
        double pL, pR, cL, cR;
        this->eos->setPressure(&pL, &rhoL, &uL);
        this->eos->setPressure(&pR, &rhoR, &uR);
        this->eos->setSoundSpeed(&cL, &rhoL, &uL);
        this->eos->setSoundSpeed(&cR, &rhoR, &uR);

        return computeHLLCFluxFromStates<dim>(rhoL, vL, uL, pL, cL,
                                              rhoR, vR, uR, pR, cR, axis);
    }
};
