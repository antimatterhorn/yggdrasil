// Copyright (C) 2025  Cody Raskin
// Kurganov–Tadmor scheme with 2nd-order MUSCL reconstruction

#pragma once
#include "gridHydroBase.hh"
#include "HLL.cc"


template<int dim>
class GridHydroKT : public GridHydroBase<dim> {
public:
    using typename GridHydroBase<dim>::Vector;
    using typename GridHydroBase<dim>::VectorField;
    using typename GridHydroBase<dim>::ScalarField;

    GridHydroKT(NodeList* nodeList, PhysicalConstants& constants,
                EquationOfState* eos, Mesh::Grid<dim>* grid)
        : GridHydroBase<dim>(nodeList, constants, eos, grid) {}

    virtual std::string name() const override { return "GridHydroKT"; }

    struct 
    ConsVars {
        double rho;
        Vector mom;
        double E;

        static ConsVars fromPrimitive(double rho, Vector v, double u) {
            ConsVars q;
            q.rho = rho;
            q.mom = v * rho;
            q.E = rho * (u + 0.5 * v.mag2());
            return q;
        }
    };

    inline double 
    minmod(double a, double b) const {
        if (a * b <= 0.0) return 0.0;
        return (std::abs(a) < std::abs(b)) ? a : b;
    }

    inline double 
    muscl_slope(double uL, double uC, double uR) const {
        return 0.5 * minmod(uC - uL, uR - uC);
    }

    inline double 
    vanleer_slope(double uL, double uC, double uR) const {
        double deltaL = uC - uL;
        double deltaR = uR - uC;

        if (deltaL * deltaR <= 0.0) return 0.0;

        return (2.0 * deltaL * deltaR) / (deltaL + deltaR);
    }

    template<typename T>
    T slopeLimitedValue(const Field<T>& field, int idx, int axis) const {
        auto neighbors = this->grid->getNeighboringCells(idx);
        int iL = neighbors[2 * axis];
        int iR = neighbors[2 * axis + 1];

        if (iL < 0 || iR < 0) return T();

        T uL = field.getValue(iL);
        T uC = field.getValue(idx);
        T uR = field.getValue(iR);

        if constexpr (std::is_same<T, double>::value) {
            return vanleer_slope(uL, uC, uR);
        } else {
            T slope;
            for (int d = 0; d < dim; ++d)
                slope[d] = vanleer_slope(uL[d], uC[d], uR[d]);
            return slope;
        }
    }

    virtual HLLFlux<dim> 
    computeFlux(int iL, int iR, int axis,
                const Field<double>& rho,
                const Field<Vector>& v,
                const Field<double>& u,
                const Field<double>& p,
                const Field<double>& cs) const override {
        double dx = this->grid->spacing(axis);

        // Cell-centered values
        double rho0L= rho.getValue(iL), rho0R = rho.getValue(iR);
        Vector v0L  = v.getValue(iL), v0R = v.getValue(iR);
        double u0L  = u.getValue(iL), u0R = u.getValue(iR);
        double p0L  = p.getValue(iL), p0R = p.getValue(iR);
        double c0L  = cs.getValue(iL), c0R = cs.getValue(iR);

        // Slopes
        double drhoL    = slopeLimitedValue(rho, iL, axis);
        double drhoR    = slopeLimitedValue(rho, iR, axis);
        Vector dvL      = slopeLimitedValue(v, iL, axis);
        Vector dvR      = slopeLimitedValue(v, iR, axis);
        double duL      = slopeLimitedValue(u, iL, axis);
        double duR      = slopeLimitedValue(u, iR, axis);
        double dpL      = slopeLimitedValue(p, iL, axis);
        double dpR      = slopeLimitedValue(p, iR, axis);
        double dcL      = slopeLimitedValue(cs, iL, axis);
        double dcR      = slopeLimitedValue(cs, iR, axis);

        // Reconstructed interface values
        double rhoL = rho0L + 0.5 * drhoL;
        double rhoR = rho0R - 0.5 * drhoR;
        Vector vL   = v0L + 0.5 * dvL;
        Vector vR   = v0R - 0.5 * dvR;
        double uL   = u0L + 0.5 * duL;
        double uR   = u0R - 0.5 * duR;
        double pL   = p0L + 0.5 * dpL;
        double pR   = p0R - 0.5 * dpR;
        double cL   = c0L + 0.5 * dcL;
        double cR   = c0R - 0.5 * dcR;

        double vnL  = vL[axis], vnR = vR[axis];

        ConsVars qL = ConsVars::fromPrimitive(rhoL, vL, uL);
        ConsVars qR = ConsVars::fromPrimitive(rhoR, vR, uR);

        HLLFlux<dim> FL, FR;
        FL.mass     = qL.rho * vnL;
        FL.momentum = qL.mom * vnL; FL.momentum[axis] += pL;
        FL.energy   = (qL.E + pL) * vnL;

        FR.mass     = qR.rho * vnR;
        FR.momentum = qR.mom * vnR; FR.momentum[axis] += pR;
        FR.energy   = (qR.E + pR) * vnR;

        double aPlus = std::max({0.0, vnL + cL, vnR + cR});
        double aMinus = std::min({0.0, vnL - cL, vnR - cR});

        double denom = aPlus - aMinus;
        HLLFlux<dim> flux;

        if (std::abs(denom) < 1e-12) {
            flux.mass       = 0.5 * (FL.mass + FR.mass);
            flux.momentum   = 0.5 * (FL.momentum + FR.momentum);
            flux.energy     = 0.5 * (FL.energy + FR.energy);
            return flux;
        }

        flux.mass       = (aPlus*FL.mass - aMinus*FR.mass + aPlus*aMinus * (qR.rho - qL.rho)) / denom;
        flux.momentum   = (FL.momentum*aPlus - FR.momentum*aMinus + aPlus*aMinus * (qR.mom - qL.mom)) / denom;
        flux.energy     = (aPlus*FL.energy - aMinus*FR.energy + aPlus*aMinus * (qR.E - qL.E)) / denom;
        return flux;
    }
};
