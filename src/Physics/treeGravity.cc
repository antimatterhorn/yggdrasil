// Copyright (C) 2026  Cody Raskin

#include "physics.hh"
#include "../Trees/spatialTree.hh"
#include <iostream>

template <int dim>
class TreeGravity : public Physics<dim> {
protected:
    double dtmin;
    double plummerLength;
    double theta;

public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    TreeGravity(NodeList* nodeList, PhysicalConstants& constants, double plummerLength, double theta = 0.5) :
        Physics<dim>(nodeList, constants),
        plummerLength(plummerLength),
        theta(theta) {

        this->template EnrollFields<Vector>({"acceleration"});
        // ACCUMULATE: contributes dvdt but does not own or finalize velocity.
        // position is read from the INTEGRATE package's sub-stage state via initialState.
        this->template EnrollStateFields<Vector>({"velocity"}, FieldPolicy::ACCUMULATE);
    }

    ~TreeGravity() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        PhysicalConstants constants = this->constants;
        int numNodes = nodeList->size();

        ScalarField* mass         = nodeList->getField<double>("mass");
        VectorField* position     = initialState->template getField<Vector>("position");
        VectorField* velocity     = initialState->template getField<Vector>("velocity");
        VectorField* acceleration = nodeList->getField<Vector>("acceleration");
        VectorField* dvdt         = deriv.template getField<Vector>("velocity");

        SpatialTree<dim> tree(position, mass);
        tree.build();

        double local_dtmin = 1e30;
        double eps2  = plummerLength;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i = 0; i < numNodes; ++i) {
            Vector a = tree.computeForceOn(i, theta, constants.G(), eps2);
            Vector v = velocity->getValue(i);

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
        double timestepCoefficient = 0.1;
        return timestepCoefficient * dtmin;
    }

    virtual std::string name() const override { return "treeGravity"; }
    virtual std::string description() const override {
        return "Barnes-Hut tree-based gravity for N-body simulations";
    }
};
