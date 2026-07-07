// Copyright (C) 2026  Cody Raskin

#include "physics.hh"
#include "../EOS/equationOfState.hh"
#include "../EOS/opacityModel.hh"
#include "../Utilities/printTable.hh"
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

class StellarEvolution : public Physics<1> {
private:
    double dm;
    double dtmin;
    int numZones;
public:
    using Vector       = Lin::Vector<1>;
    using VectorField  = Field<Vector>;
    using ScalarField  = Field<double>;

    EquationOfState* eos;
    OpacityModel* opac;
    double totalMass, centralTemperature, radius, gamma;

    StellarEvolution(NodeList* nodeList, PhysicalConstants& constants, EquationOfState* eos, OpacityModel* opac,
        double totalMass, double radius, double centralTemperature, double gamma)
        : Physics<1>(nodeList, constants), eos(eos), opac(opac), totalMass(totalMass),
        centralTemperature(centralTemperature), radius(radius), gamma(gamma),
        numZones(nodeList->size()), dtmin(1e30) {

        this->template EnrollFields<double>(
            {"pressure", "density", "mass", "radius", "specificInternalEnergy", "luminosity", "temperature"});
        this->template EnrollStateFields<double>({"specificInternalEnergy"});

        dm = totalMass / numZones;

        BuildHydrostaticModel();
    }

    void EvaluateDerivatives(const State<1>* initialState,
                             State<1>& deriv,
                             const double time,
                             const double dt) override {
        ScalarField* u     = initialState->template getField<double>("specificInternalEnergy");
        ScalarField* dUdt  = deriv.template getField<double>("specificInternalEnergy");

        NodeList* nodeList = this->nodeList;
        ScalarField* rho   = nodeList->template getField<double>("density");
        ScalarField* m     = nodeList->template getField<double>("mass");

        ScalarField T("temperature", u->size()); // local storage of temperature for computing derivatives
        eos->setTemperature(&T, rho, u);

        // Compute ε_nuc(T)
        std::vector<double> eps(u->size());
        for (int i = 0; i < u->size(); ++i)
            eps[i] = epsilonNuc(T[i]);

        // Integrate L from center
        std::vector<double> L(u->size(), 0.0);
        for (int i = 1; i < u->size(); ++i)
            L[i] = L[i-1] + eps[i-1] * dm;

        // Compute dL/dm and du/dt
        double local_dtmin = 1e30;
        for (int i = 1; i < u->size() - 1; ++i) {
            double dLdm = (L[i+1] - L[i-1]) / (m->getValue(i+1) - m->getValue(i-1));
            double du_dt = eps[i] - dLdm;
            dUdt->setValue(i, du_dt);

            if (std::abs(du_dt) > 1e-30)
                local_dtmin = std::min(local_dtmin, std::abs(u->getValue(i) / du_dt));
        }

        dUdt->setValue(0, 0.0);
        dUdt->setValue(u->size() - 1, 0.0);

        dtmin = local_dtmin;
    }

    void FinalizeStep(const State<1>* finalState) override {
        NodeList* nodeList = this->nodeList;

        ScalarField* fu = finalState->template getField<double>("specificInternalEnergy");
        ScalarField* u  = nodeList->getField<double>("specificInternalEnergy");

        u->copyValues(fu);

        #pragma omp parallel for
        for (int i = 0; i < u->size(); ++i)
            u->setValue(i, std::max((*u)[i], 1e-12));  // Prevent nonphysical negative energies

        ScalarField* rho = nodeList->getField<double>("density");
        ScalarField* T   = nodeList->getField<double>("temperature");
        ScalarField* P   = nodeList->getField<double>("pressure");

        eos->setTemperature(T, rho, u);
        eos->setPressure(P, rho, u);

        ComputeHydrostaticEquilibrium();  // Optional: maintain strict equilibrium
        ComputeLuminosity();
    }

    void ComputeHydrostaticEquilibrium() {
        NodeList* nodeList = this->nodeList;
        PhysicalConstants& constants = this->constants;

        ScalarField* P   = nodeList->getField<double>("pressure");
        ScalarField* rho = nodeList->getField<double>("density");
        ScalarField* m   = nodeList->getField<double>("mass");
        ScalarField* r   = nodeList->getField<double>("radius");

        const int nz = m->size();
        P->setValue(nz - 1, 0.0);  // Surface pressure boundary condition

        for (int i = nz - 2; i >= 0; --i) {
            double r_next = r->getValue(i+1);
            double dPdm = -constants.G() * m->getValue(i+1) / (4 * M_PI * std::pow(r_next, 4));
            double Pi = (*P)[i + 1] - dPdm * dm;
            P->setValue(i, Pi);
        }
    }

    double epsilonNuc(double T) const {
        return 1e-5 * std::pow(T / 1e7, 6);  // Toy pp-chain scaling law
    }

    void ComputeLuminosity() {
        NodeList* nodeList = this->nodeList;

        ScalarField* L = nodeList->getField<double>("luminosity");
        ScalarField* T = nodeList->getField<double>("temperature");

        (*L)[0] = 0.0;
        for (int i = 1; i < numZones; ++i) {
            double eps = epsilonNuc((*T)[i - 1]);
            (*L)[i] = (*L)[i - 1] + eps * dm;
        }
    }

    void BuildHydrostaticModel() {
        NodeList* nodeList = this->nodeList;
        PhysicalConstants& constants = this->constants;

        const int nz = numZones;
        const double tol = 1e-6;
        const int maxIter = 30;

        // Fixed spacing in radius
        const double dr = radius / nz;

        // Initial guess range for central density
        double rho_c0 = 1e4, rho_c1 = 2e4;
        double M0 = 0.0, M1 = 0.0;

        std::vector<double> best_rho, best_T, best_u, best_P, best_m, best_r, best_L;

        for (int outer = 0; outer < maxIter; ++outer) {
            double rho_c = (outer == 0 ? rho_c0 : rho_c1);

            // Init arrays
            std::vector<double> rho(nz), T(nz), u(nz), P(nz), m(nz), r(nz), L(nz, 0.0);
            r[0] = 1e-5;
            m[0] = 0.0;
            rho[0] = rho_c;
            T[0] = centralTemperature;
            eos->setInternalEnergyFromTemperature(&u[0], &rho[0], &T[0]);
            eos->setPressure(&P[0], &rho[0], &u[0]);

            for (int i = 1; i < nz; ++i) {
                r[i] = r[i-1] + dr;

                // Integrate mass
                double dm_dr = 4 * M_PI * r[i-1] * r[i-1] * rho[i-1];
                m[i] = m[i-1] + dm_dr * dr;

                // Hydrostatic equilibrium: dP/dr = -G m rho / r^2
                double dP_dr = -constants.G() * m[i-1] * rho[i-1] / (r[i-1] * r[i-1]);
                P[i] = P[i-1] + dP_dr * dr;

                // Energy generation integrated outward to this shell: dL/dr = 4 pi r^2 rho eps_nuc
                L[i] = L[i-1] + epsilonNuc(T[i-1]) * dm_dr * dr;

                // Radiative diffusion gradient: dT/dr = -L / (4 pi r^2 X), X = 4 a c T^3 / (3 kappa rho)
                double Xcond;
                opac->setConductivity(&Xcond, &rho[i-1], &T[i-1]);
                double dT_dr_rad = -L[i] / (4 * M_PI * r[i-1] * r[i-1] * Xcond);

                // Adiabatic gradient: dT/dr = (gamma-1)/gamma * (T/P) * dP/dr
                double dT_dr_ad = (gamma - 1.0) / gamma * (T[i-1] / P[i-1]) * dP_dr;

                // Schwarzschild criterion: convection sets in wherever the radiative
                // gradient needed to carry L would be steeper than the adiabatic one.
                double dT_dr = (std::abs(dT_dr_rad) > std::abs(dT_dr_ad)) ? dT_dr_ad : dT_dr_rad;

                T[i] = std::max(T[i-1] + dT_dr * dr, 1.0);
                eos->setInternalEnergyFromTemperature(&u[i], &rho[i-1], &T[i]);

                // Solve for rho_i such that P(rho_i, T[i]) = P[i]
                double rho_i = rho[i-1];
                for (int iter = 0; iter < 20; ++iter) {
                    double tmp_u, tmp_P;
                    eos->setInternalEnergyFromTemperature(&tmp_u, &rho_i, &T[i]);
                    eos->setPressure(&tmp_P, &rho_i, &tmp_u);
                    double f = tmp_P - P[i];
                    if (std::abs(f / P[i]) < tol) break;

                    // Finite difference
                    double drho = 1e-6 * rho_i;
                    double rho_pert = rho_i + drho;
                    double tmp_u_plus, tmp_P_plus;
                    eos->setInternalEnergyFromTemperature(&tmp_u_plus, &rho_pert, &T[i]);
                    eos->setPressure(&tmp_P_plus, &rho_pert, &tmp_u_plus);
                    double df = (tmp_P_plus - tmp_P) / drho;

                    rho_i -= f / df;
                    rho_i = std::max(rho_i, 1e-12);
                }

                rho[i] = rho_i;
            }

            double Mcalc = m[nz - 1];
            std::cout << "Iteration " << outer << ": rho_c = " << rho_c << ", M = " << Mcalc << "\n";

            if (std::abs((Mcalc - totalMass) / totalMass) < tol) {
                best_rho = rho;
                best_T = T;
                best_u = u;
                best_P = P;
                best_m = m;
                best_r = r;
                best_L = L;
                break;
            }

            // Update secant method
            if (outer == 0) {
                M0 = Mcalc;
                rho_c0 = rho_c;
            } else {
                M1 = Mcalc;
                double rho_new = rho_c1 - (M1 - totalMass) * (rho_c1 - rho_c0) / (M1 - M0 + 1e-12);
                rho_c0 = rho_c1;
                M0 = M1;
                rho_c1 = std::max(rho_new, 1e-4);
            }
        }

        // Rescale radius to hit target R
        double r_scale = radius / best_r.back();
        for (double& ri : best_r) ri *= r_scale;

        // Assign to fields
        ScalarField* frho = nodeList->getField<double>("density");
        ScalarField* fu   = nodeList->getField<double>("specificInternalEnergy");
        ScalarField* fP   = nodeList->getField<double>("pressure");
        ScalarField* fT   = nodeList->getField<double>("temperature");
        ScalarField* fr   = nodeList->getField<double>("radius");
        ScalarField* fm   = nodeList->getField<double>("mass");
        ScalarField* fL   = nodeList->getField<double>("luminosity");
        for (int i = 0; i < nz; ++i) {
            frho->setValue(i, best_rho[i]);
            fT->setValue(i, best_T[i]);
            fu->setValue(i, best_u[i]);
            fP->setValue(i, best_P[i]);
            fm->setValue(i, best_m[i]);
            fr->setValue(i, best_r[i]);
            fL->setValue(i, best_L[i]);
        }

        printTable(frho->size(), fr, frho, fm, fu, fP, fT, fL);
    }

    double EstimateTimestep() const override {
        double timestepCoefficient = 0.1; // fraction of the local thermal (u / du/dt) timescale
        return timestepCoefficient * dtmin;
    }

    std::string name() const override { return "StellarEvolution"; }
    std::string description() const override {
        return "1D radial stellar structure and evolution with toy pp-chain nuclear burning"; }
};
