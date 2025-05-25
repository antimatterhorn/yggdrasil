#include "physics.hh"
#include "../EOS/equationOfState.hh"
#include "../EOS/opacityModel.hh"
#include "../Mesh/grid.hh"

template <int dim>
class ThermalConduction : public Physics<dim> {
protected:
    Mesh::Grid<dim>* grid;
    EquationOfState* eos;
    OpacityModel* opac;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    ThermalConduction(NodeList* nodeList, PhysicalConstants& constants, EquationOfState* eos, OpacityModel* opac, Mesh::Grid<dim>* grid) : 
        Physics<dim>(nodeList,constants), eos(eos), grid(grid), opac(opac) {
        VerifyFields(nodeList);
    }

    virtual ~ThermalConduction() {}

    virtual void
    VerifyFields(NodeList* nodeList) {
        this->template EnrollFields<double>({"pressure", "density", "specificInternalEnergy", "soundSpeed", "temperature", "conductivity"});
        this->template EnrollStateFields<double>({"specificInternalEnergy"});
    }

    void SetConductivity() {
        NodeList* nodeList  = this->nodeList;
        int numZones        = nodeList->size();

        ScalarField* rho           = nodeList->getField<double>("density");
        ScalarField* u             = nodeList->getField<double>("specificInternalEnergy");
        ScalarField* T             = nodeList->getField<double>("temperature");
        ScalarField* x             = nodeList->getField<double>("conductivity");
        // looping and using scalar methods for speed
        for (int i = 0 ; i < numZones ; ++i) SetConductivity(rho,u,T,x,i);
    }

    void SetConductivity(ScalarField* rho, ScalarField* u, ScalarField* T, ScalarField* X, int i) {
        double rhoi = rho->getValue(i);
        double ui   = u->getValue(i);
        double Ti   = 0;
        double xi   = 0;
        eos->setTemperature(&Ti, &rhoi, &ui);
        T->setValue(i, Ti);
        opac->setConductivity(&xi, &rhoi, &Ti);
        X->setValue(i, xi);
    }

    virtual void ZeroTimeInitialize() override {
        SetConductivity();
    }

    virtual void PreStepInitialize() override {
        SetConductivity();
        this->state.updateFields(this->nodeList);
    }

    virtual void EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        int numZones = nodeList->size();

        ScalarField* rho    = nodeList->getField<double>("density");
        ScalarField* u      = initialState->template getField<double>("specificInternalEnergy");

        ScalarField* dudt   = deriv.template getField<double>("specificInternalEnergy");

        for (int i = 0 ; i < numZones ; ++i) {
            if (!grid->onBoundary(i)) {
                const Vector& ri = grid->position(i);
                double Vi = grid->cellVolume(i);
                double divFlux = 0.0;

                for (int j : grid->neighbors(i)) {
                    if (j < 0 || j >= numZones || j == i) continue;

                    const Vector& rj = grid->position(j);
                    Vector dx = rj - ri;
                    double dist2 = dx.magnitude2();
                    if (dist2 == 0.0) continue;

                    double Tij = T[j] - T[i];
                    double kij = 0.5 * (kappa[i] + kappa[j]);  // symmetric average

                    double flux = -kij * Tij / std::sqrt(dist2);  // scalar flux across face
                    double Aij = grid->faceArea(i, j);            // effective area between i and j

                    divFlux += flux * Aij;
                }

                dudt->setValue(i,divFlux / Vi);
            }
        }

        this->lastDt = dt;
    }

    virtual double EstimateTimestep() const override {
        return 1e30; // this physics package does not support setting the timestep for now
    }

    virtual void FinalizeStep(const State<dim>* finalState) override {
        PushState(finalState);
    }

    virtual void PushState(const State<dim>* stateToPush) override {
        NodeList* nodeList = this->nodeList;
        int numZones = nodeList->size();
        State<dim> state = this->state;

    }

    virtual std::string name() const override { return "ThermalConduction"; }
    virtual std::string description() const override {
        return "Thermal conduction physics"; }
};