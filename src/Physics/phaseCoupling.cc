#include "physics.hh"
#include <iostream>
#include <math>

// Kuramoto, Yoshiki (1975). H. Araki (ed.). 
// Lecture Notes in Physics, 
// International Symposium on Mathematical Problems in Theoretical Physics. 
// Vol. 39. Springer-Verlag, New York. p. 420.

template <int dim>
class SimplePhysics : public Physics<dim> {
protected:
    double couplingConstant, lightFraction;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    SimplePhysics(NodeList* nodeList,
                        PhysicalConstants& constants, 
                        double couplingConstant, double lightFraction) :
        Physics<dim>(nodeList,constants), 
        couplingConstant(couplingConstant),
        lightFraction(lightFraction) {

        this->template EnrollFields<double>({"kphase","kstrength"});
        this->template EnrollFields<Vector>({"position"});  // let kinetics handle the evo of this field
        this->template EnrollStateFields<double>({"kphase"});
    }

    ~SimplePhysics() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        auto* phase     = initialState->template getField<double>("kphase");
        auto* dph       = deriv.template getField<double>("kphase");
        auto* x         = nodeList->getField<Vector>("position");

        int numNodes = this->nodeList.size();
        for(int i=0; i<numNodes; ++i) {
            for(int j=0; j<numNodes; ++j) {
                if (j!=i) {
                    double xij = (x[i]-x[j]).magnitude();
                    double K = couplingConstant/(xij+1e-5);
                    dph[i] += K*sin(phase[j]-phase[i])/(numNodes-1);
                }
            }
        }
    }

    virtual void
    FinalChecks() override {
        NodeList* nodeList = this->nodeList;
        auto* phase     = nodeList->getField<double>("kphase");
        auto* str       = nodeList->getField<double>("kstrength");

        int numNodes = this->nodeList.size();
        for(int i=0; i<numNodes; ++i) {
            // mod phase with 2pi here
            str[i] = (fabs(sin(phase[i])) > (1.0-lightFraction));  // this is a boolean disguised as a double
        }
    }

};
