#include <complex>
#include "physics.hh"
#include "../Mesh/grid.hh"
#include <iostream>



template <int dim>
class ComplexWaveEquation : public Physics<dim> {
protected:
    Mesh::Grid<dim>* grid;
    double C;
    double dtmin;
    double dxmin = 1e30;
    std::vector<int> insideIds;
public:
    using Vector        = Lin::Vector<dim>;
    using Complex       = std::complex<double>;
    using ComplexField  = Field<Complex>;
    using ScalarField   = Field<double>;


    ComplexWaveEquation(NodeList* nodeList, PhysicalConstants& constants, Mesh::Grid<dim>* grid, double C) : 
        Physics<dim>(nodeList,constants), grid(grid), C(C) {
        VerifyWaveFields();
        grid->assignPositions(nodeList);

        auto* cs = nodeList->template getField<double>("soundSpeed");
        for (int i=0; i<nodeList->getNumNodes();++i) cs->setValue(i,C);
        for (int i=0; i<dim; ++i) dxmin = std::min(dxmin, grid->spacing(i));
    }

    void VerifyWaveFields() {
        this->template EnrollFields<Complex>({"psi"});
        this->template EnrollStateFields<Complex>({"psi"});
        this->template EnrollFields<double>({"soundSpeed", "waveEnergyDensity"});
    }

    virtual void ZeroTimeInitialize() override {
        int numNodes = this->nodeList->size();
        for (int i = 0; i < numNodes; ++i) {
            if (!grid->onBoundary(i)) insideIds.push_back(i);
        }
        this->UpdateState();
        this->InitializeBoundaries();
    }

    virtual void EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        auto* psi     = initialState->template getField<Complex>("psi");
        auto* DpsiDt  = deriv.template getField<Complex>("psi");

        auto* cs = this->nodeList->template getField<double>("soundSpeed");
        auto* e  = this->nodeList->template getField<double>("waveEnergyDensity");

        double local_dtmin = 1e30;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int p = 0; p < insideIds.size(); ++p) {
            int i = insideIds[p];

            double c = cs->getValue(i);
            Complex psi_i = psi->getValue(i);

            // Time derivative: d(psi)/dt = i * c^2 * ∇²(Re(psi)) + Im(psi)
            double phi = std::real(psi_i);
            double xi  = std::imag(psi_i);
            double laplace_phi = grid->laplacian(i, [&](int j) { return std::real(psi->getValue(j)); });

            Complex dpsi_dt(xi, c * c * laplace_phi);
            DpsiDt->setValue(i, dpsi_dt);

            // Energy density: 0.5 * (|Im(psi)|² + c² |∇Re(psi)|²)
            auto grad        = grid->gradient(i, [&](int j) { return std::real(psi->getValue(j)); });
            double grad2     = grad.mag2();

            e->setValue(i, 0.5 * (std::imag(psi_i) * std::imag(psi_i) + c * c * grad2));

            local_dtmin = std::min(local_dtmin, 0.2 * dxmin / c);
        }

        dtmin = local_dtmin;
    }

    virtual double EstimateTimestep() const override { 
        return dtmin;
    }

    virtual std::string name() const override { return "complexWaveEquation"; }
    virtual std::string description() const override {
        return "Wave equation using complex field ψ = φ + iξ";
    }

    double
    getCell(int i,int j, std::string fieldName="phi") {
        int idx = grid->index(i,j,0);
        ScalarField* phi    = this->nodeList->template getField<double>(fieldName);
        return phi->getValue(idx);
    }
};
