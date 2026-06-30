// Copyright (C) 2025  Cody Raskin

#include "kinematics.hh"
#include "../Trees/kdTree.hh"
#include <map>
#include <set>
#include <cmath>

// Discrete Element Method using the Cundall-Strack linear spring-dashpot contact model.
//
// Normal:      Fn  = max(0, kn*overlap - gammaN*vn)         (no adhesion)
// Tangential:  Ft  = -kt*deltaT - gammaT*vt, capped at mu*Fn (Coulomb)
//
// Tangential spring displacement (deltaT) is stored per contact pair and advanced
// once per timestep (timeVisited guard), so repeated RK stage calls use the same
// spring state.  In 2D the rotational surface velocity is included exactly; in 3D
// a scalar angularVelocity cannot encode the spin axis so the rolling term is omitted.

template <int dim>
class DEM : public Kinematics<dim> {
protected:
    double kn;       // normal spring stiffness     [force / length]
    double kt;       // tangential spring stiffness [force / length]
    double gammaN;   // normal damping coefficient  [force * time / length]
    double gammaT;   // tangential damping coefficient
    double mu;       // Coulomb friction coefficient (dimensionless)
    double rMax;     // cached maximum particle radius
    double timeVisited;

    using ContactKey = std::pair<int, int>;   // always stored as (i, j) with i < j
    std::map<ContactKey, Lin::Vector<dim>> tangentialDisplacement;

public:
    using Vector     = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    DEM(NodeList* nodeList, PhysicalConstants& constants,
        double kn, double kt, double gammaN, double gammaT, double mu) :
        Kinematics<dim>(nodeList, constants),
        kn(kn), kt(kt), gammaN(gammaN), gammaT(gammaT), mu(mu),
        rMax(0.0), timeVisited(-1.0) {

        this->template EnrollFields<double>({"radius", "angularVelocity", "torque"});
        this->template EnrollStateFields<double>({"angularVelocity"});
        this->template EnrollFields<Vector>({"force"});
    }

    ~DEM() {}

    virtual void
    ZeroTimeInitialize() override {
        cacheRMax();
        Kinematics<dim>::ZeroTimeInitialize();
    }

    virtual void
    PreStepInitialize() override {
        Physics<dim>::PreStepInitialize();
        cacheRMax();
    }

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        const bool isNewStep = (time != timeVisited);
        if (isNewStep) timeVisited = time;

        int n = this->nodeList->size();

        ScalarField* massFld   = this->nodeList->template getField<double>("mass");
        ScalarField* radiusFld = this->nodeList->template getField<double>("radius");
        VectorField* forceFld  = this->nodeList->template getField<Vector>("force");
        ScalarField* torqueFld = this->nodeList->template getField<double>("torque");

        VectorField* posFld   = initialState->template getField<Vector>("position");
        VectorField* velFld   = initialState->template getField<Vector>("velocity");
        ScalarField* omegaFld = initialState->template getField<double>("angularVelocity");

        // Build KDTree from current stage positions so contact search matches the state
        KDTree<dim> tree(posFld);

        for (int i = 0; i < n; ++i) {
            forceFld->setValue(i, Vector::zero());
            torqueFld->setValue(i, 0.0);
        }

        std::set<ContactKey> activeContacts;

        for (int i = 0; i < n; ++i) {
            Vector ri = posFld->getValue(i);
            double Ri = radiusFld->getValue(i);

            for (int j : tree.findNearestNeighbors(ri, Ri + rMax)) {
                if (j <= i) continue;   // process each unordered pair once

                Vector rj      = posFld->getValue(j);
                double Rj      = radiusFld->getValue(j);
                Vector rij     = ri - rj;            // points from j to i
                double dist    = rij.magnitude();
                double overlap = Ri + Rj - dist;
                if (overlap <= 0.0) continue;

                Vector nHat  = rij / dist;           // unit normal, j → i
                Vector vi    = velFld->getValue(i);
                Vector vj    = velFld->getValue(j);
                Vector vRel  = vi - vj;
                double vRelN = vRel * nHat;

                // --- Normal spring-dashpot (no tensile contact) ---
                double Fn     = std::max(0.0, kn * overlap - gammaN * vRelN);
                Vector Fn_vec = Fn * nHat;

                // --- Tangential relative velocity at contact surface ---
                Vector vRelTang = vRel - vRelN * nHat;
                if constexpr (dim == 2) {
                    // Exact rotational contribution: ω × r_contact projected onto tangent
                    // Surface velocity of i at contact: -Ri * omegaI * t̂
                    // Surface velocity of j at contact: +Rj * omegaJ * t̂
                    // Relative tangential slip adds: (-Ri*omegaI - Rj*omegaJ) * t̂
                    Vector tHat  = {-nHat[1], nHat[0]};   // CCW perpendicular to nHat
                    double vRoll = -Ri * omegaFld->getValue(i) - Rj * omegaFld->getValue(j);
                    vRelTang     = vRelTang + vRoll * tHat;
                }

                // --- Cundall-Strack tangential spring ---
                ContactKey key = {i, j};
                if (isNewStep) activeContacts.insert(key);

                // try_emplace: inserts zero vector only if key is new; leaves existing unchanged
                auto result  = tangentialDisplacement.try_emplace(key, Vector::zero());
                Vector& deltaT = result.first->second;

                if (isNewStep) {
                    // Re-project onto current contact plane (corrects for rotating contact)
                    deltaT = deltaT - (deltaT * nHat) * nHat;
                    // Advance by tangential slip
                    deltaT = deltaT + vRelTang * dt;
                }

                Vector Ft    = -kt * deltaT - gammaT * vRelTang;
                double FtMag = Ft.magnitude();

                // Coulomb limit: cap force and correct stored spring displacement
                if (FtMag > mu * Fn && FtMag > 0.0) {
                    Ft = Ft * (mu * Fn / FtMag);
                    if (isNewStep)
                        deltaT = (-1.0 / kt) * (Ft + gammaT * vRelTang);
                }

                // --- Accumulate forces (Newton's 3rd law) ---
                forceFld->setValue(i, forceFld->getValue(i) + Fn_vec + Ft);
                forceFld->setValue(j, forceFld->getValue(j) - Fn_vec - Ft);

                // --- Torque (scalar z-component of r × F_tangential) ---
                // τ_i = -Ri * (n̂ × F_t)_z;  same formula for j (both spin same way)
                // Only physically meaningful in 2D; in 3D scalar ω has no defined axis.
                if constexpr (dim == 2) {
                    double crossNFt = nHat[0] * Ft[1] - nHat[1] * Ft[0];
                    torqueFld->setValue(i, torqueFld->getValue(i) - Ri * crossNFt);
                    torqueFld->setValue(j, torqueFld->getValue(j) - Rj * crossNFt);
                }
            }
        }

        // Prune contacts that are no longer active (once per timestep)
        if (isNewStep) {
            for (auto it = tangentialDisplacement.begin(); it != tangentialDisplacement.end(); ) {
                if (activeContacts.find(it->first) == activeContacts.end())
                    it = tangentialDisplacement.erase(it);
                else
                    ++it;
            }
        }

        // --- Fill time derivatives ---
        VectorField* dxdt = deriv.template getField<Vector>("position");
        VectorField* dvdt = deriv.template getField<Vector>("velocity");
        ScalarField* dwdt = deriv.template getField<double>("angularVelocity");

        const double Ifac = (dim == 2) ? 0.5 : 0.4;   // disk (2D) or sphere (3D)

        for (int i = 0; i < n; ++i) {
            double mi   = massFld->getValue(i);
            double Ri   = radiusFld->getValue(i);
            double Ii   = Ifac * mi * Ri * Ri;

            dxdt->setValue(i, velFld->getValue(i));
            dvdt->setValue(i, forceFld->getValue(i) / mi);
            dwdt->setValue(i, torqueFld->getValue(i) / Ii);
        }

        this->lastDt = dt;
    }

    virtual double
    EstimateTimestep() const override {
        ScalarField* massFld = this->nodeList->template getField<double>("mass");
        int n = this->nodeList->size();
        double mMin = 1e30;
        for (int i = 0; i < n; ++i)
            mMin = std::min(mMin, massFld->getValue(i));
        // 1/10 of the half-period of the lightest spring-mass oscillator
        return 0.1 * M_PI * std::sqrt(mMin / kn);
    }

    virtual std::string name()        const override { return "DEM"; }
    virtual std::string description() const override {
        return "Discrete Element Method with Cundall-Strack linear spring-dashpot contact model";
    }

private:
    void cacheRMax() {
        ScalarField* radiusFld = this->nodeList->template getField<double>("radius");
        int n = this->nodeList->size();
        rMax = 0.0;
        for (int i = 0; i < n; ++i)
            rMax = std::max(rMax, radiusFld->getValue(i));
    }
};
