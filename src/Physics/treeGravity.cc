// Copyright (C) 2025  Cody Raskin

#include "kinematics.hh"
#include "../Trees/spatialTree.hh"
#include <iostream>

template <int dim>
class TreeGravity : public Kinematics<dim> {
protected:
    double dtmin;
    double plummerLength;

public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;



    TreeGravity(NodeList* nodeList, PhysicalConstants& constants, double plummerLength) :
        Kinematics<dim>(nodeList, constants),
        plummerLength(plummerLength) {}

    ~TreeGravity() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        PhysicalConstants constants = this->constants;
        int numNodes = nodeList->size();

        ScalarField* mass           = nodeList->getField<double>("mass");
        VectorField* position       = initialState->template getField<Vector>("position");
        VectorField* acceleration   = nodeList->getField<Vector>("acceleration");
        VectorField* velocity       = initialState->template getField<Vector>("velocity");
        VectorField* dxdt           = deriv.template getField<Vector>("position");
        VectorField* dvdt           = deriv.template getField<Vector>("velocity");

        // Build tree
        SpatialTree<dim> tree(position, mass);
        tree.build();

        double local_dtmin = 1e30;
        double theta = 0.5;
        double eps2 = plummerLength;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i = 0; i < numNodes; ++i) {
            Vector a = tree.computeForceOn(i, theta, constants.G(), eps2);
            Vector v = velocity->getValue(i);

            acceleration->setValue(i, a);
            dxdt->setValue(i, v + dt * a);
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
        double timestepCoefficient = 1e-2;
        return timestepCoefficient * std::sqrt(dtmin);
    }

    virtual std::string name() const override { return "treeGravity"; }
    virtual std::string description() const override {
        return "Barnes-Hut tree-based gravity for N-body simulations";
    }
};