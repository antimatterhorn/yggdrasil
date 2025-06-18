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
    EvaluateDerivatives(const State<2>* initialState, State<2>& deriv, const double time, const double dt) override {
        auto* c1 = initialState->template getField<double>("c1");
        auto* c2 = initialState->template getField<double>("c2");
        auto* c3 = initialState->template getField<double>("c3");

        auto* dc1 = deriv.template getField<double>("c1");
        auto* dc2 = deriv.template getField<double>("c2");
        auto* dc3 = deriv.template getField<double>("c3");

        auto* rho = deriv.template getField<double>("density");

        const int sx = grid->size_x();
        const int sy = grid->size_y();

        #pragma omp parallel for collapse(2)
        for (int j = 0; j < sy; ++j) {
            for (int i = 0; i < sx; ++i) {
                int idx = grid->index(i, j);

                double u[3] = { c1->getValue(idx), c2->getValue(idx), c3->getValue(idx) };
                double dudt[3] = { 0.0, 0.0, 0.0 };

                // Compute total density
                double r = u[0] + u[1] + u[2];
                rho->setValue(idx, r);

                // Diffusion term for each component
                for (int c = 0; c < 3; ++c) {
                    double del = 0.0;
                    for (int dj = -1; dj <= 1; ++dj) {
                        for (int di = -1; di <= 1; ++di) {
                            int ii = (i + di + sx) % sx;
                            int jj = (j + dj + sy) % sy;
                            int nidx = grid->index(ii, jj);

                            double fac = (di != 0 && dj != 0) ? 0.05 : 0.2;
                            del += fac * (
                                (c == 0 ? c1->getValue(nidx) :
                                (c == 1 ? c2->getValue(nidx) : c3->getValue(nidx)))
                            );
                        }
                    }

                    // Logistic and cyclic reaction terms
                    double competitor = u[(c + 1) % 3];
                    dudt[c] = D * del + u[c] * (1.0 - r) - A * u[c] * competitor;
                }

                // Store result in derivative state
                dc1->setValue(idx, dudt[0]);
                dc2->setValue(idx, dudt[1]);
                dc3->setValue(idx, dudt[2]);
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
