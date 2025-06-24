// Copyright (C) 2025  Cody Raskin

#include "integrator.hh"

template <int dim>
class CrankNicolsonIntegrator : public Integrator<dim> {
public:
    CrankNicolsonIntegrator(std::vector<Physics<dim>*> packages, double dtmin, bool verbose = false) :
        Integrator<dim>(packages, dtmin, verbose) {}

    ~CrankNicolsonIntegrator() {}

    virtual State<dim> 
    Integrate(Physics<dim>* physics) override { 
        double dt = this->dt;
        double time = this->time;

        const State<dim>* state = physics->getState();

        State<dim> k1(state->size());
        State<dim> k2(state->size());

        k1.ghost(state);
        k2.ghost(state);

        // Derivatives at current time/state
        physics->EvaluateDerivatives(state, k1, time, 0.0);

        // Initial guess: Euler forward
        State<dim> predicted = state->deepCopy();
        predicted += k1 * dt;

        // Fixed-point iteration
        const double tolerance = 1e-10;
        const int maxIterations = 10;

        for (int iter = 0; iter < maxIterations; ++iter) {
            physics->EvaluateDerivatives(&predicted, k2, time, dt);

            State<dim> newPredicted = state->deepCopy();
            State<dim> avgDeriv = k1.deepCopy();
            avgDeriv += k2;
            avgDeriv *= 0.5 * dt;

            newPredicted += avgDeriv;

            double delta = (newPredicted - predicted).L2Norm();

            if (this->verbose)
                std::cout << "CrankNicolson iteration " << iter << ": Δ = " << delta << "\n";

            if (delta < tolerance) {
                this->dtMultiplier *= (iter < maxIterations*0.5 ? 
                        1.2 : (iter > maxIterations*0.8 ? 0.8 : 1));
                break;
            }                    

            predicted.swap(newPredicted);  // update for next iteration
        }

        return predicted;
    }
};
