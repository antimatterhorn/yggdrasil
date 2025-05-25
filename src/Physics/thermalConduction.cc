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
        grid->assignPositions(nodeList);
    }

    virtual ~ThermalConduction() {}

    virtual void
    VerifyFields(NodeList* nodeList) {
        this->template EnrollFields<double>({"pressure", "density", "specificInternalEnergy", "soundSpeed", "temperature", "conductivity"});
        this->template EnrollFields<Vector>({"position"});
        this->template EnrollStateFields<double>({"specificInternalEnergy"});
    }

    void SetConductivity() {
        NodeList* nodeList  = this->nodeList;
        int numZones        = nodeList->size();

        ScalarField* rho           = nodeList->getField<double>("density");
        ScalarField* u             = nodeList->getField<double>("specificInternalEnergy");
        ScalarField* T             = nodeList->getField<double>("temperature");
        ScalarField* X             = nodeList->getField<double>("conductivity");
        // looping and using scalar methods for speed
        for (int i = 0 ; i < numZones ; ++i) SetConductivity(rho,u,T,X,i);
    }

    void SetConductivity(ScalarField* rho, ScalarField* u, ScalarField* T, ScalarField* X, int i) {
        double rhoi = rho->getValue(i);
        double ui   = u->getValue(i);
        double Ti   = T->getValue(i);
        double Xi   = X->getValue(i);
        eos->setTemperature(&Ti, &rhoi, &ui);
        T->setValue(i, Ti);
        opac->setConductivity(&Xi, &rhoi, &Ti);
        X->setValue(i, Xi);
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
        ScalarField* T      = nodeList->getField<double>("temperature");
        ScalarField* X      = nodeList->getField<double>("conductivity");

        ScalarField* dudt   = deriv.template getField<double>("specificInternalEnergy");

        double local_dtmin = 1e30;
        double dx2 = grid->getdx() * grid->getdx();  // assume uniform dx for now

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i = 0 ; i < numZones ; ++i) {
            if (!grid->onBoundary(i)) {
                const Vector& ri = grid->getPosition(i);
                double Vi = grid->cellVolume(i);
                double divFlux = 0.0;
                double Xi = X->getValue(i);
                double Ti = T->getValue(i);

                SetConductivity(rho,u,T,X,i);

                for (int j : grid->neighbors(i)) {
                    if (j < 0 || j >= numZones || j == i) continue;

                    double Xj = X->getValue(j);
                    double Tj = T->getValue(j);

                    const Vector& rj = grid->getPosition(j);
                    Vector dx = rj - ri;
                    double dist2 = dx.mag2();
                    if (dist2 == 0.0) continue;

                    double Tij = Tj-Ti;
                    double Xij = 0.5 * (Xi+Xj);  // symmetric average

                    double flux = -Xij * Tij / std::sqrt(dist2);  // scalar flux across face
                    double Aij = grid->faceArea(i, j);            // effective area between i and j

                    divFlux += flux * Aij;
                }

                dudt->setValue(i,divFlux / Vi);

                double rhoi = rho->getValue(i);  // density

                // Approximate cv = du/dT numerically:
                double ui = u->getValue(i);
                double dT = Ti * 1.0e-4 + 1e-10;  // small perturbation
                double ui_plus, Ti_plus = Ti + dT;
                eos->setInternalEnergy(&ui_plus, &rhoi, &Ti_plus);
                double cv = (ui_plus - ui) / dT;

                double D = Xi / (rhoi * cv + 1e-30);  // prevent divide by zero
                double dt_candidate = 0.5 * dx2 / (D + 1e-30);

                local_dtmin = std::min(local_dtmin, dt_candidate);
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

        ScalarField* u      = nodeList->getField<double>("specificInternalEnergy");
        ScalarField* fu     = stateToPush->template getField<double>("specificInternalEnergy");

        u->copyValues(fu);

        if (stateToPush != &(state))
        {
            ScalarField* su = state.template getField<double>("specificInternalEnergy");
            su->copyValues(fu);
        }
    }

    double getCell(int i, int j, const std::string& fieldName = "pressure") const {
        int idx = grid->index(i, j, 0);
        ScalarField* field = this->nodeList->template getField<double>(fieldName);
        return field->getValue(idx);
    }

    virtual std::string name() const override { return "ThermalConduction"; }
    virtual std::string description() const override {
        return "Thermal conduction physics"; }
};