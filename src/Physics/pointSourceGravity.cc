// Copyright (C) 2026  Cody Raskin

#include "physics.hh"
#include <iostream>

template <int dim>
class PointSourceGravity : public Physics<dim> {
protected:
    Lin::Vector<dim> pointSourceLocation;
    Lin::Vector<dim> pointSourceVelocity;
    double pointSourceMass;
    double dtmin;
    double lastUpdateTime = -1.0;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    PointSourceGravity(NodeList* nodeList,
                        PhysicalConstants& constants,
                        Vector& pointSourceLocation,
                        Vector& pointSourceVelocity,
                        double pointSourceMass) :
        Physics<dim>(nodeList,constants),
        pointSourceLocation(pointSourceLocation),
        pointSourceVelocity(pointSourceVelocity),
        pointSourceMass(pointSourceMass) {

        this->template EnrollFields<Vector>({"acceleration"});
        // ACCUMULATE: contributes dvdt but does not own or finalize velocity.
        // position is read from the INTEGRATE package's sub-stage state via initialState.
        this->template EnrollStateFields<Vector>({"velocity"}, FieldPolicy::ACCUMULATE);
    }

    ~PointSourceGravity() {}

    Vector getPointSourceLocation() const { return pointSourceLocation; }
    void setPointSourceLocation(Vector loc) { pointSourceLocation = loc; }

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        // Advance point source once per timestep (not per RK sub-stage).
        if (time != lastUpdateTime) {
            pointSourceLocation += pointSourceVelocity * dt;
            lastUpdateTime = time;
        }

        NodeList* nodeList = this->nodeList;
        PhysicalConstants constants = this->constants;
        int numNodes = nodeList->size();

        VectorField* position     = initialState->template getField<Vector>("position");
        VectorField* velocity     = initialState->template getField<Vector>("velocity");
        VectorField* acceleration = nodeList->getField<Vector>("acceleration");
        VectorField* dvdt         = deriv.template getField<Vector>("velocity");

        double local_dtmin = 1e30;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i=0; i<numNodes; ++i) {
            Vector pos = position->getValue(i);
            Vector r   = pointSourceLocation - pos;
            Vector a   = pointSourceMass * constants.G() / r.mag2() * r.normal();
            Vector v   = velocity->getValue(i);
            acceleration->setValue(i, a);
            dvdt->setValue(i, a);

            double amag = a.mag2();
            double vmag = v.mag2();
            if (amag > 0.0)
                local_dtmin = std::min(local_dtmin, vmag / amag);
        }

        dtmin = local_dtmin;
        this->lastDt = dt;
    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 1e-4;
        return timestepCoefficient * dtmin;
    }

    virtual std::string name() const override { return "pointSourceGravity"; }
    virtual std::string description() const override {
        return "Point source gravity physics package for particles"; }
};
