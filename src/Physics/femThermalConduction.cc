// Copyright (C) 2025  Cody Raskin

#include "physics.hh"
#include "../Mesh/femesh.hh"
#include "../EOS/opacityModel.hh"

template <int dim>
class FEMHeatConduction : public Physics<dim> {

private:
    Mesh::FEMesh<dim>* mesh;
    OpacityModel* opac;
    double dtmin;

public:
    using Vector = Lin::Vector<dim>;
    using ScalarField = Field<double>;

    FEMHeatConduction(NodeList* nodeList, PhysicalConstants& constants,
                      Mesh::FEMesh<dim>* mesh, OpacityModel* opac)
        : Physics<dim>(nodeList, constants), mesh(mesh), opac(opac) {
        this->template EnrollFields<double>({"temperature", "density", "conductivity"});
        this->template EnrollStateFields<double>({"temperature"});
    }

    void EvaluateDerivatives(const State<dim>* state, State<dim>& deriv,
                             const double time, const double dt) override {
        const auto& elements     = mesh->getElements();
        const auto& connectivity = mesh->getConnectivityMap();
        const auto& positions    = mesh->getNodes();

        ScalarField* temperature  = state->template getField<double>("temperature");
        ScalarField* dTdt         = deriv.template getField<double>("temperature");

        ScalarField* density      = this->nodeList->template getField<double>("density");
        ScalarField* conductivity = this->nodeList->template getField<double>("conductivity");

        opac->setConductivity(conductivity, density, temperature);

        int numNodes = this->nodeList->size();
        std::vector<double> rhs(numNodes, 0.0);
        std::vector<double> mass(numNodes, 0.0);

        for (size_t e = 0; e < elements.size(); ++e) {
            const auto& elem = *elements[e];
            const auto& conn = connectivity[e];

            // Element-averaged conductivity
            double k = 0.0;
            for (size_t i = 0; i < conn.size(); ++i)
                k += conductivity->getValue(conn[i]);
            k /= static_cast<double>(conn.size());

            Eigen::MatrixXd Ke = k * elem.computeStiffnessMatrix(positions);
            Eigen::VectorXd Me = elem.computeLumpedMassMatrix(positions);

            Eigen::VectorXd Te(conn.size());
            for (size_t i = 0; i < conn.size(); ++i)
                Te(i) = temperature->getValue(conn[i]);

            Eigen::VectorXd fe = -Ke * Te;

            for (size_t i = 0; i < conn.size(); ++i) {
                rhs[conn[i]]  += fe(i);
                mass[conn[i]] += Me(i);
            }
        }

        double local_dtmin = 1e30;
        for (int i = 0; i < numNodes; ++i) {
            double dT = (mass[i] > 0.0) ? rhs[i] / mass[i] : 0.0;
            dTdt->setValue(i, dT);
            double k = conductivity->getValue(i);
            if (k > 0.0 && mass[i] > 0.0)
                local_dtmin = std::min(local_dtmin, 0.5 * mass[i] / k);
        }
        dtmin = local_dtmin;
    }

    virtual double
    EstimateTimestep() const override { return dtmin; }

    std::string name() const override { return "FEMHeatConduction"; }
    std::string description() const override { return "FEM heat conduction with opacity-based conductivity"; }
};
