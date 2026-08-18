// Copyright (C) 2026  Cody Raskin

#include "physics.hh"
#include "../Trees/fmmTree.hh"
#include <iostream>

template <int dim>
class FMMGravity : public Physics<dim> {
protected:
    double dtmin;
    double plummerLength;
    int maxSourcesPerLeaf;

public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    FMMGravity(NodeList* nodeList, PhysicalConstants& constants, double plummerLength, int maxSourcesPerLeaf = 16) :
        Physics<dim>(nodeList, constants),
        plummerLength(plummerLength),
        maxSourcesPerLeaf(maxSourcesPerLeaf) {

        this->template EnrollFields<Vector>({"acceleration"});
        this->template EnrollStateFields<Vector>({"velocity"}, FieldPolicy::ACCUMULATE);
    }

    ~FMMGravity() {}

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

        FMMTree<dim> tree(position, mass, maxSourcesPerLeaf);
        tree.build();
        tree.upwardPass();
        tree.downwardPass(constants.G());

        double local_dtmin = 1e30;
        double eps2 = plummerLength;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i = 0; i < numNodes; ++i) {
            Vector a = tree.computeForceOn(i, constants.G(), eps2);
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

    virtual std::string name() const override { return "fmmGravity"; }
    virtual std::string description() const override {
        return "Fast multipole method (quadrupole) gravity for N-body simulations";
    }
};
