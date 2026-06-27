// Copyright (C) 2025  Cody Raskin

#include "physics.hh"
#include "../Mesh/femesh.hh"
#include "../Materials/material.hh"
#include <cmath>

template <int dim>
class FEMLinearElasticity : public Physics<dim> {
protected:
    Mesh::FEMesh<dim>* mesh;
    FEMMaterial<dim>*  material;
    double rho;
    Lin::Vector<dim>   bodyForce;
    double alpha, beta;  // Rayleigh damping: f_damp = -(alpha*M + beta*K_diag) * v
    double dtmin = 1e30;
    std::vector<Lin::Vector<dim>> refPositions;

public:
    using Vector      = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    FEMLinearElasticity(NodeList* nodeList, PhysicalConstants& constants,
                        Mesh::FEMesh<dim>* mesh, FEMMaterial<dim>* material,
                        double rho, Vector bodyForce,
                        double alpha = 0.0, double beta = 0.0)
        : Physics<dim>(nodeList, constants),
          mesh(mesh), material(material),
          rho(rho), bodyForce(bodyForce),
          alpha(alpha), beta(beta) {
        this->template EnrollFields<Vector>({"position", "velocity", "displacement"});
        this->template EnrollFields<double>({"vonMises", "sigmaXX", "sigmaYY", "sigmaXY"});
        this->template EnrollStateFields<Vector>({"position", "velocity"});
    }

    void ZeroTimeInitialize() override {
        const auto& meshNodes = mesh->getNodes();
        refPositions.resize(meshNodes.size());
        VectorField* posField = this->nodeList->template getField<Vector>("position");
        for (size_t i = 0; i < meshNodes.size(); ++i) {
            refPositions[i] = meshNodes[i];
            posField->setValue(i, meshNodes[i]);
        }

        assembleNodalMass();
        this->UpdateState();
        this->InitializeBoundaries();
    }

    void EvaluateDerivatives(const State<dim>* state, State<dim>& deriv,
                             const double time, const double dt) override {
        const auto& elements     = mesh->getElements();
        const auto& connectivity = mesh->getConnectivityMap();
        const auto& meshPos      = mesh->getNodes();

        VectorField* pos  = state->template getField<Vector>("position");
        VectorField* vel  = state->template getField<Vector>("velocity");
        VectorField* dxdt = deriv.template getField<Vector>("position");
        VectorField* dvdt = deriv.template getField<Vector>("velocity");
        ScalarField* massField = this->nodeList->template getField<double>("mass");

        Eigen::MatrixXd D = material->materialMatrix();
        int numNodes = this->nodeList->size();

        std::vector<double> rhs(numNodes * dim, 0.0);
        std::vector<double> kdiag(numNodes, 0.0);

        for (size_t e = 0; e < elements.size(); ++e) {
            const auto& conn = connectivity[e];
            int ndofs = static_cast<int>(conn.size()) * dim;

            Eigen::MatrixXd Ke = elements[e]->computeStructuralStiffnessMatrix(meshPos, D);

            Eigen::VectorXd ue(ndofs);
            for (size_t i = 0; i < conn.size(); ++i) {
                Vector u = pos->getValue(conn[i]) - refPositions[conn[i]];
                for (int d = 0; d < dim; ++d)
                    ue(i * dim + d) = u[d];
            }

            Eigen::VectorXd fe = -Ke * ue;

            for (size_t i = 0; i < conn.size(); ++i) {
                for (int d = 0; d < dim; ++d)
                    rhs[conn[i] * dim + d] += fe(i * dim + d);
                for (int d = 0; d < dim; ++d)
                    kdiag[conn[i]] += Ke(i * dim + d, i * dim + d);
            }
        }

        double local_dtmin = 1e30;

        for (int i = 0; i < numNodes; ++i) {
            double mi = massField->getValue(i);
            if (mi <= 0.0) continue;

            Vector v = vel->getValue(i);

            Vector f_elastic;
            for (int d = 0; d < dim; ++d)
                f_elastic[d] = rhs[i * dim + d];

            double c = alpha * mi + beta * kdiag[i];
            Vector f_damping;
            for (int d = 0; d < dim; ++d)
                f_damping[d] = -c * v[d];

            Vector a = (f_elastic + f_damping) * (1.0 / mi) + bodyForce;

            dxdt->setValue(i, v + dt * a);
            dvdt->setValue(i, a);

            if (kdiag[i] > 0.0)
                local_dtmin = std::min(local_dtmin, std::sqrt(mi / kdiag[i]));
        }

        dtmin = local_dtmin;
        this->lastDt = dt;
    }

    void FinalChecks() override {
        const auto& elements     = mesh->getElements();
        const auto& connectivity = mesh->getConnectivityMap();
        const auto& meshPos      = mesh->getNodes();

        VectorField* pos      = this->nodeList->template getField<Vector>("position");
        VectorField* dispField = this->nodeList->template getField<Vector>("displacement");
        int numNodes = this->nodeList->size();

        for (int i = 0; i < numNodes; ++i)
            dispField->setValue(i, pos->getValue(i) - refPositions[i]);

        if constexpr (dim == 2) {
            Eigen::MatrixXd D = material->materialMatrix();

            ScalarField* sxxField = this->nodeList->template getField<double>("sigmaXX");
            ScalarField* syyField = this->nodeList->template getField<double>("sigmaYY");
            ScalarField* sxyField = this->nodeList->template getField<double>("sigmaXY");
            ScalarField* vmField  = this->nodeList->template getField<double>("vonMises");

            std::vector<double> sxx(numNodes, 0.0), syy(numNodes, 0.0), sxy(numNodes, 0.0);
            std::vector<int> count(numNodes, 0);

            for (size_t e = 0; e < elements.size(); ++e) {
                const auto& conn = connectivity[e];
                int ndofs = static_cast<int>(conn.size()) * dim;

                Eigen::VectorXd ue(ndofs);
                for (size_t i = 0; i < conn.size(); ++i) {
                    Vector u = pos->getValue(conn[i]) - refPositions[conn[i]];
                    for (int d = 0; d < dim; ++d)
                        ue(i * dim + d) = u[d];
                }

                Eigen::VectorXd strain = elements[e]->computeStrain(meshPos, ue);
                Eigen::VectorXd sigma  = D * strain;

                for (size_t i = 0; i < conn.size(); ++i) {
                    sxx[conn[i]] += sigma(0);
                    syy[conn[i]] += sigma(1);
                    sxy[conn[i]] += sigma(2);
                    count[conn[i]]++;
                }
            }

            for (int i = 0; i < numNodes; ++i) {
                if (count[i] == 0) continue;
                double sx = sxx[i] / count[i];
                double sy = syy[i] / count[i];
                double txy = sxy[i] / count[i];
                sxxField->setValue(i, sx);
                syyField->setValue(i, sy);
                sxyField->setValue(i, txy);
                vmField->setValue(i, std::sqrt(sx*sx - sx*sy + sy*sy + 3.0*txy*txy));
            }
        }
    }

    double EstimateTimestep() const override { return 0.5 * dtmin; }

    std::string name()        const override { return "FEMLinearElasticity"; }
    std::string description() const override { return "Linear elastic FEM structural mechanics"; }

private:
    void assembleNodalMass() {
        const auto& elements     = mesh->getElements();
        const auto& connectivity = mesh->getConnectivityMap();
        const auto& meshPos      = mesh->getNodes();

        ScalarField* massField = this->nodeList->template getField<double>("mass");
        int numNodes = this->nodeList->size();
        for (int i = 0; i < numNodes; ++i) massField->setValue(i, 0.0);

        for (size_t e = 0; e < elements.size(); ++e) {
            const auto& conn = connectivity[e];
            Eigen::VectorXd Me = rho * elements[e]->computeLumpedMassMatrix(meshPos);
            for (size_t i = 0; i < conn.size(); ++i)
                massField->setValue(conn[i], massField->getValue(conn[i]) + Me(i));
        }
    }
};
