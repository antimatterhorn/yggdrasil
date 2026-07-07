// Copyright (C) 2026  Cody Raskin

#include "integrator.hh"

template <int dim>
class RungeKutta2Integrator : public Integrator<dim> {
protected:

public:
    RungeKutta2Integrator(std::vector<Physics<dim>*> packages, double dtmin, bool verbose = false) :
        Integrator<dim>(packages,dtmin,verbose) {}

    ~RungeKutta2Integrator() {}

    virtual State<dim>
    Integrate(Physics<dim>* physics,
              const std::vector<Physics<dim>*>& accumulators) override {
        double dt = this->dt;
        double time = this->time;

        const State<dim>* state = physics->getState();
        State<dim> interim = state->deepCopy();
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

        evalWithAccum(state, k1, time, 0);

        interim += k1 * dt;

        evalWithAccum(&interim, k2, time, dt);

        State<dim> newState = state->deepCopy();

        k1 += k2;
        k1 *= 0.5 * dt;
        newState += k1;

        return newState;
    }
};
