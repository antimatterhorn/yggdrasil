// Copyright (C) 2026  Cody Raskin

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
    double dtmin;
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
        this->template EnrollStateFields<double>({"specificInternalEnergy"}, FieldPolicy::INTEGRATE);
    }

    void SetConductivity() {
        SetConductivity(this->nodeList->template getField<double>("specificInternalEnergy"));
    }

    // Refresh temperature and conductivity from the given energy field, which may
    // belong to an RK sub-stage State rather than to the NodeList.
    void SetConductivity(ScalarField* u) {
        int numZones = this->nodeList->size();

        ScalarField* rho           = this->nodeList->template getField<double>("density");
        ScalarField* T             = this->nodeList->template getField<double>("temperature");
        ScalarField* X             = this->nodeList->template getField<double>("conductivity");
        // looping and using scalar methods for speed
        #pragma omp parallel for
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
        this->UpdateState();
        this->InitializeBoundaries();
    }

    virtual void PreStepInitialize() override {
        SetConductivity();
        this->state.updateFields(this->nodeList);
    }

    // EvaluateDerivatives leaves T and X at the last RK sub-stage; resync them
    // against the finalized energy so anything reading the NodeList sees the step's
    // actual end-of-step values.
    virtual void FinalChecks() override {
        SetConductivity();
    }

    virtual void EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        int numZones = this->nodeList->size();

        ScalarField* rho    = this->nodeList->template getField<double>("density");
        ScalarField* u      = initialState->template getField<double>("specificInternalEnergy");
        ScalarField* T      = this->nodeList->template getField<double>("temperature");
        ScalarField* X      = this->nodeList->template getField<double>("conductivity");

        ScalarField* dudt   = deriv.template getField<double>("specificInternalEnergy");

        double local_dtmin = 1e30;
        double dx2 = grid->getdx() * grid->getdx();  // assume uniform dx for now

        SetConductivity(u);

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i = 0 ; i < numZones ; ++i) {
            if (!grid->onBoundary(i)) {
                const Vector& ri = grid->getPosition(i);
                double Vi = grid->cellVolume(i);
                double divFlux = 0.0;
                double Xi = X->getValue(i);
                double Ti = T->getValue(i);

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

                    divFlux -= flux * Aij;
                }

                double rhoi = rho->getValue(i);  // density

                dudt->setValue(i, divFlux / (rhoi * Vi));

                // Approximate cv = du/dT numerically:
                double ui = u->getValue(i);
                double dT = std::max(std::abs(Ti) * 1e-4, 1e-10);
                double ui_plus, Ti_plus = Ti + dT;
                eos->setInternalEnergyFromTemperature(&ui_plus, &rhoi, &Ti_plus);
                double cv = (ui_plus - ui) / dT;

                // Clamp or skip invalid results
                if (cv <= 0.0 || std::isnan(cv) || std::isinf(cv)) continue;

                double D = Xi / (rhoi * cv);
                if (D <= 0.0 || std::isnan(D) || std::isinf(D)) continue;

                double dt_candidate = 0.5 * dx2 / D;
                if (dt_candidate > 0.0)
                    local_dtmin = std::min(local_dtmin, dt_candidate);
            }
        }
        dtmin = local_dtmin;

        //std::cout << "dtmin = " << dtmin << std::endl;
        this->lastDt = dt;
    }

    virtual double EstimateTimestep() const override {
        double timestepCoefficient = 0.25; // Adjust as needed
        double timestep = timestepCoefficient * dtmin;

        return timestep;
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