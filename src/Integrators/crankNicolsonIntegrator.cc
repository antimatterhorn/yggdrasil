// Copyright (C) 2026  Cody Raskin

#include "integrator.hh"

template <int dim>
class CrankNicolsonIntegrator : public Integrator<dim> {
public:
    CrankNicolsonIntegrator(std::vector<Physics<dim>*> packages, double dtmin, bool verbose = false) :
        Integrator<dim>(packages, dtmin, verbose) {}

    ~CrankNicolsonIntegrator() {}

    virtual State<dim>
    Integrate(Physics<dim>* physics,
              const std::vector<Physics<dim>*>& accumulators) override {
        double dt = this->dt;
        double time = this->time;

        const State<dim>* state = physics->getState();

        State<dim> k1(state->size());
        State<dim> k2(state->size());

        k1.ghost(state);
        k2.ghost(state);

        auto evalWithAccum = [&](const State<dim>* s, State<dim>& k, double t, double subDt) {
            physics->EvaluateDerivatives(s, k, t, subDt);
            for (auto* acc : accumulators) {
                State<dim> accK(state->size());
                accK.ghost(acc->getState());
                acc->EvaluateDerivatives(s, accK, t, subDt);
                k.accumulateFrom(accK, acc->accumulateFieldNames());
            }
        };

        evalWithAccum(state, k1, time, 0.0);

        // Initial guess: Euler forward. Both sides of the convergence test below
        // must have their rims refreshed alike, or the rim contributes a constant
        // residual and the fixed-point iteration never converges.
        State<dim> predicted = state->deepCopy();
        predicted += k1 * dt;
        physics->ApplyStageBoundaries(&predicted);

        // Fixed-point iteration
        const double tolerance = 1e-10;
        const int maxIterations = 10;

        for (int iter = 0; iter < maxIterations; ++iter) {
            evalWithAccum(&predicted, k2, time, dt);

            State<dim> newPredicted = state->deepCopy();
            State<dim> avgDeriv = k1.deepCopy();
            avgDeriv += k2;
            avgDeriv *= 0.5 * dt;

            newPredicted += avgDeriv;
            physics->ApplyStageBoundaries(&newPredicted);

            double delta = (newPredicted - predicted).L2Norm();

            if (this->verbose)
                std::cout << "CrankNicolson iteration " << iter << ": Δ = " << delta << "\n";

            if (delta < tolerance) {
                this->dtMultiplier *= (iter < maxIterations * 0.5 ?
                        1.2 : (iter > maxIterations * 0.8 ? 0.8 : 1));
                break;
            }

            predicted.swap(newPredicted);
        }

        return predicted;
    }
};
