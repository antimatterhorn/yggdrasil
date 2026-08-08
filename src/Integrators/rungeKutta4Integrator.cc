// Copyright (C) 2026  Cody Raskin

#include "integrator.hh"

template <int dim>
class RungeKutta4Integrator : public Integrator<dim> {
protected:

public:
    RungeKutta4Integrator(std::vector<Physics<dim>*> packages, double dtmin, bool verbose = false) :
        Integrator<dim>(packages,dtmin,verbose) {}

    ~RungeKutta4Integrator() {}

    virtual State<dim>
    Integrate(Physics<dim>* physics,
              const std::vector<Physics<dim>*>& accumulators) override {
        double dt = this->dt;
        double time = this->time;

        const State<dim>* state = physics->getState();
        State<dim> interim = state->deepCopy();
        State<dim> k1(state->size());
        State<dim> k2(state->size());
        State<dim> k3(state->size());
        State<dim> k4(state->size());

        k1.ghost(state);
        k2.ghost(state);
        k3.ghost(state);
        k4.ghost(state);

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

        interim += k1 * (dt / 2.0);
        physics->ApplyStageBoundaries(&interim);
        evalWithAccum(&interim, k2, time, dt / 2.0);

        interim.clone(state);
        interim += k2 * (dt / 2.0);
        physics->ApplyStageBoundaries(&interim);
        evalWithAccum(&interim, k3, time, dt / 2.0);

        interim.clone(state);
        interim += k3 * dt;
        physics->ApplyStageBoundaries(&interim);
        evalWithAccum(&interim, k4, time, dt);

        k2 *= 2.0;
        k3 *= 2.0;
        k1 += k2;
        k1 += k3;
        k1 += k4;
        k1 *= (dt / 6.0);

        State<dim> newState = state->deepCopy();
        newState += k1;
        return newState;
    }
};
