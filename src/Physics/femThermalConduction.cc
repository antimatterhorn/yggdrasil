#include "physics.hh"
#include "../Mesh/femesh.hh"

template <int dim>
class FEMHeatConduction : public Physics<dim> {

private:
    FEMesh<dim>* mesh;

public:
    using Vector = Lin::Vector<dim>;
    using ScalarField = Field<double>;

    FEMHeatConduction(NodeList* nodeList, PhysicalConstants& constants, Mesh::FEMesh<dim>* mesh)
        : Physics<dim>(nodeList,constants), mesh(mesh) {
        this->template EnrollFields<double>({"temperature","mass"});
        this->themplate EnrollStateFields<double>({"temperature"});
    }

    void EvaluateDerivatives(const State<dim>* state, State<dim>& deriv, const double time, const double dt) override {
        const auto& elements        = mesh->getElements();
        const auto& connectivity    = mesh->getConnectivityMap();
        const auto& positions       = mesh->getNodes();

        ScalarField* temperature = state->template getField<double>("temperature");
        ScalarField* dTdt        = deriv.template getField<double>("temperature");

        ScalarField rhs("rhs", this->nodeList);
        ScalarField mass("mass", this->nodeList);
        rhs.setAll(0.0);
        mass.setAll(0.0);

        for (size_t e = 0; e < elements.size(); ++e) {
            const auto& elem = *elements[e];
            const auto& conn = connectivity[e];

            Eigen::MatrixXd Ke = elem.computeStiffnessMatrix(positions);
            Eigen::VectorXd Me = elem.computeLumpedMassMatrix(positions);

            Eigen::VectorXd Te(conn.size());
            for (size_t i = 0; i < conn.size(); ++i) {
                Te(i) = temperature->getValue(conn[i]);
            }

            Eigen::VectorXd fe = -Ke * Te;

            for (size_t i = 0; i < conn.size(); ++i) {
                rhs[conn[i]] += fe(i);
                mass[conn[i]] += Me(i);
            }
        }

        for (size_t i = 0; i < dTdt->size(); ++i) {
            dTdt->setValue(i, rhs[i] / mass[i]);
        }
    }

    std::string name() const override { return "FEMHeatConduction"; }

};
