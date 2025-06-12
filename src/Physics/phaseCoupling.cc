#include "physics.hh"
#include <iostream>
#include <cmath>
#include "../Trees/kdTree.hh"

// Kuramoto, Yoshiki (1975). H. Araki (ed.). 
// Lecture Notes in Physics, 
// International Symposium on Mathematical Problems in Theoretical Physics. 
// Vol. 39. Springer-Verlag, New York. p. 420.

template <int dim>
class PhaseCoupling : public Physics<dim> {
protected:
    double couplingConstant, lightFraction;
    double searchRadius=0;
    inline double mod2pi(double x) {
        double twopi = 2.0 * M_PI;
        double result = std::fmod(x, twopi);
        return result < 0 ? result + twopi : result;
    }
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    PhaseCoupling(NodeList* nodeList,
                        PhysicalConstants& constants, 
                        double couplingConstant, double lightFraction) :
        Physics<dim>(nodeList,constants), 
        couplingConstant(couplingConstant),
        lightFraction(lightFraction) {
        Enroll();
    }

    PhaseCoupling(NodeList* nodeList,
                        PhysicalConstants& constants, 
                        double couplingConstant, double lightFraction,
                        double searchRadius) :
        Physics<dim>(nodeList,constants), 
        couplingConstant(couplingConstant),
        lightFraction(lightFraction),
        searchRadius(searchRadius) {
        Enroll();
    }

    void
    Enroll() {
        this->template EnrollFields<double>({"kphase","kstrength","komega"});
        this->template EnrollFields<Vector>({"position"});  // let kinetics handle the evo of this field
        this->template EnrollStateFields<double>({"kphase"});
    }

    ~PhaseCoupling() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        auto* phase     = initialState->template getField<double>("kphase");
        auto* dph       = deriv.template getField<double>("kphase");
        auto* x         = nodeList->getField<Vector>("position");
        auto* omega     = nodeList->getField<double>("komega");

        int numNodes = nodeList->size();
        KDTree tree(x);  // this isn't zero cost, but it also isn't extremely expensive
        #pragma omp parallel for
        for(int i=0; i<numNodes; ++i) {
            Vector xi = x->getValue(i);
            dph->setValue(i, omega->getValue(i));
            std::vector<int> neighbors  = tree.findNearestNeighbors(xi, searchRadius);
            int nbrs = (searchRadius == 0? numNodes : neighbors.size());
            double norm = (searchRadius == 0? (numNodes-1) : std::max(1, (int)neighbors.size()));
            for(int j=0; j<nbrs; ++j) {
                int k = (searchRadius == 0? j : neighbors[j]);
                double pd = pdot(phase, x, i, k);
                dph->setValue(i, dph->getValue(i) + pd/norm);
            }
        }
    }

    double pdot(ScalarField *phase, VectorField *pos, int i, int j) {
        Vector xi = pos->getValue(i);
        Vector xj = pos->getValue(j);
        double xij = std::max((xi - xj).magnitude(), 0.005);
        double K = couplingConstant / xij;
        double pi = phase->getValue(i);
        double pj = phase->getValue(j);

        return K*sin(pj-pi);
    }

    virtual void
    FinalChecks() override {
        NodeList* nodeList = this->nodeList;
        auto* phase     = nodeList->getField<double>("kphase");
        auto* str       = nodeList->getField<double>("kstrength");

        int numNodes = nodeList->size();
        for(int i=0; i<numNodes; ++i) {
            phase->setValue(i, mod2pi(phase->getValue(i)));
            double strength = (lightFraction == 0 ? sin(phase->getValue(i)) : (fabs(sin(phase->getValue(i))) > (1.0-lightFraction)));
            str->setValue(i,strength);  // this is a boolean disguised as a double
        }
    }

};
