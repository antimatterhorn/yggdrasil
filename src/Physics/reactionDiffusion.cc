// Copyright (C) 2025  Cody Raskin

#include "physics.hh"
#include "../Mesh/grid.hh"
#include <iostream>
#include <cmath>

// A rock-paper-scissors chemical diffusion model

class ReactionDiffusion : public Physics<2> {
protected:
    Mesh::Grid<2>* grid;
    double A,D;
    std::vector<int> insideIds;
public:
    using Vector = Lin::Vector<2>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    ReactionDiffusion(NodeList* nodeList,
                        PhysicalConstants& constants, Mesh::Grid<2>* grid, 
                        double A, double D) :
        Physics<2>(nodeList,constants), A(A), D(D), grid(grid) {
        Enroll();
        grid->assignPositions(nodeList);
    }

    void
    Enroll() {
        this->template EnrollFields<double>({"density","c1","c2","c3"});
        this->template EnrollFields<Vector>({"position"});
        this->template EnrollStateFields<double>({"density","c1","c2","c3"});
    }

    ~ReactionDiffusion() {}

    virtual void
    ZeroTimeInitialize() override {
        int numNodes = this->nodeList->size();
        for (int i=0; i<numNodes; ++i) {
            if (!grid->onBoundary(i))
                insideIds.push_back(i);
        }

        this->UpdateState();
        this->InitializeBoundaries();
    }

    virtual void
    EvaluateDerivatives(const State<2>* initialState, State<2>& deriv, const double time, const double dt) override {
        auto* c1 = initialState->template getField<double>("c1");
        auto* c2 = initialState->template getField<double>("c2");
        auto* c3 = initialState->template getField<double>("c3");

        auto* dc1 = deriv.template getField<double>("c1");
        auto* dc2 = deriv.template getField<double>("c2");
        auto* dc3 = deriv.template getField<double>("c3");

        auto* rho = deriv.template getField<double>("density");

        const int nx = grid->size_x();
        const int ny = grid->size_y();
        const double dx = grid->dx;
        const double dy = grid->dy;

        double* c1_data = c1->data();
        double* c2_data = c2->data();
        double* c3_data = c3->data();

        double* dc1_data = dc1->data();
        double* dc2_data = dc2->data();
        double* dc3_data = dc3->data();

        double* rho_data = rho->data();
        const int total = c1->size();  // total number of grid cells

        const int N = insideIds.size();
        const int* insideIds_data = insideIds.data();

        #pragma omp target data map(to: insideIds_data[0:N], c1_data[0:total], c2_data[0:total], c3_data[0:total]) \
                        map(from: dc1_data[0:total], dc2_data[0:total], dc3_data[0:total], rho_data[0:total])
        {
            #pragma omp target teams distribute parallel for
            for (int p = 0; p < N; ++p) {
                int idx = insideIds_data[p];
                int ix = idx % nx;
                int iy = idx / nx;

                double u[3] = {
                    c1_data[idx],
                    c2_data[idx],
                    c3_data[idx]
                };
                double dudt[3] = {0.0, 0.0, 0.0};
                double r = u[0] + u[1] + u[2];
                rho_data[idx] = r;

                for (int c = 0; c < 3; ++c) {
                    double del = 0.0;
                    for (int dj = -1; dj <= 1; ++dj) {
                        for (int di = -1; di <= 1; ++di) {
                            if (di == 0 && dj == 0) continue;
                            double fac = (di != 0 && dj != 0) ? 0.05 : 0.2;
                            int nidx = (iy + dj) * nx + (ix + di);
                            del += fac * (c == 0 ? c1_data[nidx] :
                                        c == 1 ? c2_data[nidx] :
                                                c3_data[nidx]);
                        }
                    }
                    dudt[c] = D * del + u[c] * (1.0 - r) - A * u[c] * u[(c + 1) % 3];
                }

                dc1_data[idx] = dudt[0];
                dc2_data[idx] = dudt[1];
                dc3_data[idx] = dudt[2];
            }
        }

    }

    virtual double
    EstimateTimestep() const override {
        double dt_reaction = 1.0 / (A + 1e-6);
        double dt_diffusion = 1.0 / (4.0 * D + 1e-6);  // rough bound for stability
        return std::min(dt_reaction, dt_diffusion);
    }

    virtual void
    FinalChecks() override {
        NodeList* nodeList = this->nodeList;
        auto* c1 = nodeList->getField<double>("c1");
        auto* c2 = nodeList->getField<double>("c2");
        auto* c3 = nodeList->getField<double>("c3");

        for (int i = 0; i < nodeList->size(); ++i) {
            double r = c1->getValue(i) + c2->getValue(i) + c3->getValue(i);
            if (r > 1.1 || r < -0.1) {
                c1->setValue(i,c1->getValue(i)/r);
                c2->setValue(i,c2->getValue(i)/r);
                c3->setValue(i,c3->getValue(i)/r);
                //std::cerr << "Warning: total concentration out of bounds at node " << i << ": " << r << "\n";
            }
        }
    }

    std::array<double, 3>
    getCell(int i, int j, std::string fieldName = "doesntmatter") {
        NodeList* nodeList = this->nodeList;
        auto* c1 = nodeList->getField<double>("c1");
        auto* c2 = nodeList->getField<double>("c2");
        auto* c3 = nodeList->getField<double>("c3");
        int idx  = grid->index(i,j);
        return std::array<double, 3>{
            c1->getValue(idx),
            c2->getValue(idx),
            c3->getValue(idx)
        };
    }

    virtual std::string name() const override { return "reactionDiffusion"; }
    virtual std::string description() const override {
        return "A 3-body reaction diffusion model"; }

};
