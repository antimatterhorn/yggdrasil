// Copyright (C) 2026  Cody Raskin

#pragma once
#include "physics.hh"

template<int dim>
struct HLLFlux {
    double mass;
    Lin::Vector<dim> momentum;
    double energy;
};

template<int dim>
HLLFlux<dim>
computeHLLEFlux(int iL, int iR, int axis,
               const Field<double>& rho,
               const Field<Lin::Vector<dim>>& v,
               const Field<double>& u,
               const Field<double>& p,
               const Field<double>& cs) {
    using Vector = Lin::Vector<dim>;

    // Left state
    double rhoL = rho.getValue(iL);
    Vector vL = v.getValue(iL);
    double uL = u.getValue(iL);
    double pL = p.getValue(iL);
    double cL = cs.getValue(iL);
    Vector momL = vL * rhoL;
    double eL = uL + 0.5 * vL.mag2();
    double EL = rhoL * eL;

    // Right state
    double rhoR = rho.getValue(iR);
    Vector vR = v.getValue(iR);
    double uR = u.getValue(iR);
    double pR = p.getValue(iR);
    double cR = cs.getValue(iR);
    Vector momR = vR * rhoR;
    double eR = uR + 0.5 * vR.mag2();
    double ER = rhoR * eR;

    // Normal velocities
    double vnL = vL[axis];
    double vnR = vR[axis];

    // Fluxes for left and right states
    double massFluxL = rhoL * vnL;
    Vector momFluxL = momL * vnL;
    momFluxL[axis] += pL;
    double energyFluxL = vnL * (EL + pL);

    double massFluxR = rhoR * vnR;
    Vector momFluxR = momR * vnR;
    momFluxR[axis] += pR;
    double energyFluxR = vnR * (ER + pR);

    // Estimate wave speeds
    double sL = std::min(vnL - cL, vnR - cR);
    double sR = std::max(vnL + cL, vnR + cR);

    // HLL flux
    HLLFlux<dim> result;

    if (sL >= 0.0) {
        result.mass = massFluxL;
        result.momentum = momFluxL;
        result.energy = energyFluxL;
    } else if (sR <= 0.0) {
        result.mass = massFluxR;
        result.momentum = momFluxR;
        result.energy = energyFluxR;
    } else {
        double sRSum = sR - sL;
        result.mass     = (sR * massFluxL - sL * massFluxR + sR * sL * (rhoR - rhoL)) / sRSum;
        result.momentum = (sR * momFluxL - sL * momFluxR + sR * sL * (momR - momL)) / sRSum;
        result.energy   = (sR * energyFluxL - sL * energyFluxR + sR * sL * (ER - EL)) / sRSum;
    }

    return result;
}

// Normal-vector, state-based form of the same two-wave HLL flux, for solvers
// with no grid axis to key off of (ALEMeshHydroHLLE). Unlike HLLC's
// star-region solve, HLL has a single averaged middle state, so there's no
// analogue of the denomL/denomR sign bug that computeHLLCFluxFromStates had.
template<int dim>
HLLFlux<dim>
computeHLLEFluxFromStates(
    double rhoL, const Lin::Vector<dim>& vL, double uL, double pL, double cL,
    double rhoR, const Lin::Vector<dim>& vR, double uR, double pR, double cR,
    const Lin::Vector<dim>& n) {

    using Vector = Lin::Vector<dim>;

    Vector momL = rhoL * vL;
    Vector momR = rhoR * vR;
    double eL = uL + 0.5 * vL.mag2();
    double eR = uR + 0.5 * vR.mag2();
    double EL = rhoL * eL;
    double ER = rhoR * eR;

    double vnL = vL.dot(n);
    double vnR = vR.dot(n);

    double massFluxL = rhoL * vnL;
    Vector momFluxL = momL * vnL + pL * n;
    double energyFluxL = vnL * (EL + pL);

    double massFluxR = rhoR * vnR;
    Vector momFluxR = momR * vnR + pR * n;
    double energyFluxR = vnR * (ER + pR);

    double sL = std::min(vnL - cL, vnR - cR);
    double sR = std::max(vnL + cL, vnR + cR);

    HLLFlux<dim> result;

    if (sL >= 0.0) {
        result.mass = massFluxL;
        result.momentum = momFluxL;
        result.energy = energyFluxL;
    } else if (sR <= 0.0) {
        result.mass = massFluxR;
        result.momentum = momFluxR;
        result.energy = energyFluxR;
    } else {
        double sRSum = sR - sL;
        result.mass     = (sR * massFluxL - sL * massFluxR + sR * sL * (rhoR - rhoL)) / sRSum;
        result.momentum = (sR * momFluxL - sL * momFluxR + sR * sL * (momR - momL)) / sRSum;
        result.energy   = (sR * energyFluxL - sL * energyFluxR + sR * sL * (ER - EL)) / sRSum;
    }

    return result;
}

// ALE (moving-face) form: face moves with normal speed S = w.n. The flux an
// observer riding the face sees is F(U) - S*U, but only the *advective* part
// of each flux component shifts by S -- the pressure terms (p*n in momentum,
// p*vn in energy) use the true (unshifted) vn, since pressure-work doesn't
// transform the same way advection does. Region membership is decided by
// comparing the face's own speed S against the lab-frame wave speeds sL/sR
// (rather than against 0, as in the static case), since we want to know
// which part of the Riemann fan the moving face itself sits in. The single
// HLL-averaged conserved state (rho_hll etc) is a lab-frame Rankine-Hugoniot
// consequence of the two outer waves and does not depend on S; the flux
// correction term uses the face-relative wave speed (sL - S). At S=0 this
// reduces algebraically to exactly the static formula above.
template<int dim>
HLLFlux<dim>
computeHLLEFluxFromStatesALE(
    double rhoL, const Lin::Vector<dim>& vL, double uL, double pL, double cL,
    double rhoR, const Lin::Vector<dim>& vR, double uR, double pR, double cR,
    const Lin::Vector<dim>& n, double S) {

    using Vector = Lin::Vector<dim>;

    Vector momL = rhoL * vL;
    Vector momR = rhoR * vR;
    double eL = uL + 0.5 * vL.mag2();
    double eR = uR + 0.5 * vR.mag2();
    double EL = rhoL * eL;
    double ER = rhoR * eR;

    double vnL = vL.dot(n);
    double vnR = vR.dot(n);

    double sL = std::min(vnL - cL, vnR - cR);
    double sR = std::max(vnL + cL, vnR + cR);

    // Lab-frame (unshifted) static fluxes, needed only to build the
    // HLL-averaged intermediate state below.
    double massFluxL = rhoL * vnL;
    Vector momFluxL  = momL * vnL + pL * n;
    double energyFluxL = vnL * (EL + pL);

    double massFluxR = rhoR * vnR;
    Vector momFluxR  = momR * vnR + pR * n;
    double energyFluxR = vnR * (ER + pR);

    double massFluxL_ALE = rhoL * (vnL - S);
    Vector momFluxL_ALE  = momL * (vnL - S) + pL * n;
    double energyFluxL_ALE = EL * (vnL - S) + pL * vnL;

    double massFluxR_ALE = rhoR * (vnR - S);
    Vector momFluxR_ALE  = momR * (vnR - S) + pR * n;
    double energyFluxR_ALE = ER * (vnR - S) + pR * vnR;

    HLLFlux<dim> result;

    if (S <= sL) {
        result.mass = massFluxL_ALE;
        result.momentum = momFluxL_ALE;
        result.energy = energyFluxL_ALE;
    } else if (S >= sR) {
        result.mass = massFluxR_ALE;
        result.momentum = momFluxR_ALE;
        result.energy = energyFluxR_ALE;
    } else {
        double sRSum = sR - sL;
        double rho_hll = (sR * rhoR - sL * rhoL - massFluxR + massFluxL) / sRSum;
        Vector mom_hll = (sR * momR - sL * momL - momFluxR + momFluxL) / sRSum;
        double E_hll   = (sR * ER   - sL * EL   - energyFluxR + energyFluxL) / sRSum;

        result.mass     = massFluxL_ALE     + (sL - S) * (rho_hll - rhoL);
        result.momentum = momFluxL_ALE      + (sL - S) * (mom_hll - momL);
        result.energy   = energyFluxL_ALE   + (sL - S) * (E_hll - EL);
    }

    return result;
}

//LEGACY UNLIMITED VERSION
template<int dim>
HLLFlux<dim>
computeHLLCFlux(int iL, int iR, int axis,
                const Field<double>& rho,
                const Field<Lin::Vector<dim>>& v,
                const Field<double>& u,
                const Field<double>& p,
                const Field<double>& cs) {
    using Vector = Lin::Vector<dim>;

    constexpr double tiny = 1e-12;
    constexpr double rho_floor = 1e-12;

    // Left state
    double rhoL = std::max(rho.getValue(iL), rho_floor);
    Vector vL   = v.getValue(iL);
    double uL   = std::max(u.getValue(iL), rho_floor);
    double pL   = p.getValue(iL);
    double cL   = cs.getValue(iL);
    Vector momL = rhoL * vL;
    double eL   = uL + 0.5 * vL.mag2();
    double EL   = rhoL * eL;
    double vnL  = vL[axis];

    // Right state
    double rhoR = std::max(rho.getValue(iR), rho_floor);
    Vector vR   = v.getValue(iR);
    double uR   = std::max(u.getValue(iR), rho_floor);
    double pR   = p.getValue(iR);
    double cR   = cs.getValue(iR);
    Vector momR = rhoR * vR;
    double eR   = uR + 0.5 * vR.mag2();
    double ER   = rhoR * eR;
    double vnR  = vR[axis];

    // Estimate wave speeds
    double sL = std::min(vnL - cL, vnR - cR);
    double sR = std::max(vnL + cL, vnR + cR);

    // Estimate contact wave speed
    double numerator   = pR - pL + rhoL * vnL * (sL - vnL) - rhoR * vnR * (sR - vnR);
    double denominator = rhoL * (sL - vnL) - rhoR * (sR - vnR);
    if (std::abs(denominator) < tiny)
        denominator = (denominator >= 0.0) ? tiny : -tiny;
    double sStar = numerator / denominator;

    // Compute fluxes
    double massFluxL = rhoL * vnL;
    Vector momFluxL = vnL * momL;
    momFluxL[axis] += pL;
    double energyFluxL = vnL * (EL + pL);

    double massFluxR = rhoR * vnR;
    Vector momFluxR = vnR * momR;
    momFluxR[axis] += pR;
    double energyFluxR = vnR * (ER + pR);

    HLLFlux<dim> flux;

    if (sL >= 0.0) {
        flux.mass     = massFluxL;
        flux.momentum = momFluxL;
        flux.energy   = energyFluxL;
        return flux;
    }

    if (sR <= 0.0) {
        flux.mass     = massFluxR;
        flux.momentum = momFluxR;
        flux.energy   = energyFluxR;
        return flux;
    }

    // Star region calculations
    double denomL = std::max(sL - sStar, tiny);
    double denomR = std::max(sR - sStar, tiny);

    double rhoSL = std::max(rhoL * (sL - vnL) / denomL, rho_floor);
    double rhoSR = std::max(rhoR * (sR - vnR) / denomR, rho_floor);

    Vector n = unitAxis<dim>(axis);
    Vector vL_star = vL - n * vnL;
    Vector vR_star = vR - n * vnR;

    Vector momSL = momL + (sStar - vnL) * (rhoSL * vL_star);
    Vector momSR = momR + (sStar - vnR) * (rhoSR * vR_star);

    double ESL = EL + (sStar - vnL) * (rhoSL * (eL + 0.5 * (sStar - vnL)));
    double ESR = ER + (sStar - vnR) * (rhoSR * (eR + 0.5 * (sStar - vnR)));

    if (sStar >= 0.0) {
        flux.mass     = massFluxL + sL * (rhoSL - rhoL);
        flux.momentum = momFluxL + sL * (momSL - momL);
        flux.energy   = energyFluxL + sL * (ESL - EL);
    } else {
        flux.mass     = massFluxR + sR * (rhoSR - rhoR);
        flux.momentum = momFluxR + sR * (momSR - momR);
        flux.energy   = energyFluxR + sR * (ESR - ER);
    }

    // Final check: if any flux component is NaN, print warning
    if (!std::isfinite(flux.mass) ||
        !std::isfinite(flux.energy) ||
        !std::isfinite(flux.momentum[0])) {
        std::cerr << "HLLC NaN: sStar=" << sStar
                  << ", rhoL=" << rhoL << ", rhoR=" << rhoR
                  << ", vnL=" << vnL << ", vnR=" << vnR
                  << ", sL=" << sL << ", sR=" << sR
                  << ", denom=" << denominator << std::endl;
    }

    return flux;
}

//LIMITED VERSION
// Normal-vector form: the actual Riemann solve, independent of any grid axis.
// Reused by GridHydroHLLC (via the axis-based overload below, which just
// passes unitAxis<dim>(axis)) and by ALEMeshHydroHLLC (which has a genuine
// face normal, not an axis, and calls this directly).
template<int dim>
HLLFlux<dim>
computeHLLCFluxFromStates(
    double rhoL, const Lin::Vector<dim>& vL, double uL, double pL, double cL,
    double rhoR, const Lin::Vector<dim>& vR, double uR, double pR, double cR,
    const Lin::Vector<dim>& n) {

    using Vector = Lin::Vector<dim>;

    Vector momL = rhoL * vL;
    Vector momR = rhoR * vR;
    double eL = uL + 0.5 * vL.mag2();
    double eR = uR + 0.5 * vR.mag2();
    double EL = rhoL * eL;
    double ER = rhoR * eR;
    double vnL = vL.dot(n);
    double vnR = vR.dot(n);

    // Estimate wave speeds
    double sL = std::min(vnL - cL, vnR - cR);
    double sR = std::max(vnL + cL, vnR + cR);

    // Contact wave speed estimate (Toro)
    double numerator = pR - pL + rhoL * vnL * (sL - vnL) - rhoR * vnR * (sR - vnR);
    double denominator = rhoL * (sL - vnL) - rhoR * (sR - vnR);
    const double tiny = 1e-12;
    if (std::abs(denominator) < tiny)
        denominator = (denominator >= 0.0) ? tiny : -tiny;
    double sStar = numerator / denominator;

    // Left flux
    double massFluxL = rhoL * vnL;
    Vector momFluxL = vnL * momL + pL * n;
    double energyFluxL = vnL * (EL + pL);

    // Right flux
    double massFluxR = rhoR * vnR;
    Vector momFluxR = vnR * momR + pR * n;
    double energyFluxR = vnR * (ER + pR);

    HLLFlux<dim> flux;
    if (sL >= 0.0) {
        flux.mass = massFluxL;
        flux.momentum = momFluxL;
        flux.energy = energyFluxL;
    } else if (sR <= 0.0) {
        flux.mass = massFluxR;
        flux.momentum = momFluxR;
        flux.energy = energyFluxR;
    } else {
        // Toro's HLLC star states U*_K, with F*_K = F_K + s_K (U*_K - U_K). The
        // star density factor rho_K (s_K - vn_K)/(s_K - sStar) is positive:
        // numerator and denominator share a sign (both negative for K=L, both
        // positive for K=R), so the s_K - sStar guard is magnitude-based to keep
        // that sign. (n is already a parameter here, not built from an axis.)
        double dsL = sL - sStar;
        double dsR = sR - sStar;
        if (std::abs(dsL) < tiny) dsL = (dsL >= 0.0) ? tiny : -tiny;
        if (std::abs(dsR) < tiny) dsR = (dsR >= 0.0) ? tiny : -tiny;

        double rhoSL = rhoL * (sL - vnL) / dsL;
        double rhoSR = rhoR * (sR - vnR) / dsR;

        // Normal velocity becomes sStar, tangential unchanged.
        Vector momSL = rhoSL * (vL + n * (sStar - vnL));
        Vector momSR = rhoSR * (vR + n * (sStar - vnR));

        // Total energy density, Toro eq. 10.39 (eK = EK/rhoK).
        double ESL = rhoSL * (eL + (sStar - vnL) * (sStar + pL / (rhoL * (sL - vnL))));
        double ESR = rhoSR * (eR + (sStar - vnR) * (sStar + pR / (rhoR * (sR - vnR))));

        if (sStar >= 0.0) {
            flux.mass = massFluxL + sL * (rhoSL - rhoL);
            flux.momentum = momFluxL + sL * (momSL - momL);
            flux.energy = energyFluxL + sL * (ESL - EL);
        } else {
            flux.mass = massFluxR + sR * (rhoSR - rhoR);
            flux.momentum = momFluxR + sR * (momSR - momR);
            flux.energy = energyFluxR + sR * (ESR - ER);
        }
    }

    return flux;
}

// ALE (moving-face) form of HLLC, face moving with normal speed S = w.n.
// Same principle as computeHLLEFluxFromStatesALE above: sL, sR, sStar and
// the star-region conserved states are lab-frame quantities (Rankine-Hugoniot
// across each wave, unaffected by how fast the face moves), computed exactly
// as in the static solve. What changes is (a) the base left/right flux, which
// picks up the ALE correction F(U)-S*U (advective terms shift by S, pressure
// terms don't -- naively shifting vn inside the p*vn energy term too would
// introduce a spurious -S*p error), and (b) which lab-frame speed the
// face's own trajectory S is compared against to pick a region, and the wave
// speed used in the star-region flux correction term, both shifted by -S. At
// S=0 this reduces algebraically to exactly the static formula above.
template<int dim>
HLLFlux<dim>
computeHLLCFluxFromStatesALE(
    double rhoL, const Lin::Vector<dim>& vL, double uL, double pL, double cL,
    double rhoR, const Lin::Vector<dim>& vR, double uR, double pR, double cR,
    const Lin::Vector<dim>& n, double S) {

    using Vector = Lin::Vector<dim>;

    Vector momL = rhoL * vL;
    Vector momR = rhoR * vR;
    double eL = uL + 0.5 * vL.mag2();
    double eR = uR + 0.5 * vR.mag2();
    double EL = rhoL * eL;
    double ER = rhoR * eR;
    double vnL = vL.dot(n);
    double vnR = vR.dot(n);

    double sL = std::min(vnL - cL, vnR - cR);
    double sR = std::max(vnL + cL, vnR + cR);

    double numerator = pR - pL + rhoL * vnL * (sL - vnL) - rhoR * vnR * (sR - vnR);
    double denominator = rhoL * (sL - vnL) - rhoR * (sR - vnR);
    const double tiny = 1e-12;
    if (std::abs(denominator) < tiny)
        denominator = (denominator >= 0.0) ? tiny : -tiny;
    double sStar = numerator / denominator;

    double massFluxL_ALE = rhoL * (vnL - S);
    Vector momFluxL_ALE  = momL * (vnL - S) + pL * n;
    double energyFluxL_ALE = EL * (vnL - S) + pL * vnL;

    double massFluxR_ALE = rhoR * (vnR - S);
    Vector momFluxR_ALE  = momR * (vnR - S) + pR * n;
    double energyFluxR_ALE = ER * (vnR - S) + pR * vnR;

    HLLFlux<dim> flux;
    if (S <= sL) {
        flux.mass = massFluxL_ALE;
        flux.momentum = momFluxL_ALE;
        flux.energy = energyFluxL_ALE;
    } else if (S >= sR) {
        flux.mass = massFluxR_ALE;
        flux.momentum = momFluxR_ALE;
        flux.energy = energyFluxR_ALE;
    } else {
        double denomL = std::min(sL - sStar, -tiny);
        double denomR = std::max(sR - sStar, tiny);

        double rhoSL = std::max(rhoL * (sL - vnL) / denomL, tiny);
        double rhoSR = std::max(rhoR * (sR - vnR) / denomR, tiny);

        Vector vL_star = vL - n * vnL;
        Vector vR_star = vR - n * vnR;

        Vector momSL = momL + (sStar - vnL) * (rhoSL * vL_star);
        Vector momSR = momR + (sStar - vnR) * (rhoSR * vR_star);

        double ESL = EL + (sStar - vnL) * (rhoSL * (eL + 0.5 * (sStar - vnL)));
        double ESR = ER + (sStar - vnR) * (rhoSR * (eR + 0.5 * (sStar - vnR)));

        if (S <= sStar) {
            flux.mass     = massFluxL_ALE   + (sL - S) * (rhoSL - rhoL);
            flux.momentum = momFluxL_ALE    + (sL - S) * (momSL - momL);
            flux.energy   = energyFluxL_ALE + (sL - S) * (ESL - EL);
        } else {
            flux.mass     = massFluxR_ALE   + (sR - S) * (rhoSR - rhoR);
            flux.momentum = momFluxR_ALE    + (sR - S) * (momSR - momR);
            flux.energy   = energyFluxR_ALE + (sR - S) * (ESR - ER);
        }
    }

    return flux;
}

// Axis-aligned form, kept for GridHydroHLLC: identical to the general form
// above with n = unitAxis<dim>(axis).
template<int dim>
HLLFlux<dim>
computeHLLCFluxFromStates(
    double rhoL, const Lin::Vector<dim>& vL, double uL, double pL, double cL,
    double rhoR, const Lin::Vector<dim>& vR, double uR, double pR, double cR,
    int axis) {
    return computeHLLCFluxFromStates<dim>(rhoL, vL, uL, pL, cL,
                                          rhoR, vR, uR, pR, cR,
                                          unitAxis<dim>(axis));
}
