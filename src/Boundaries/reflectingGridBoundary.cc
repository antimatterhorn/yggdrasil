// Copyright (C) 2026  Cody Raskin

#include <vector>
#include <algorithm>
#include <unordered_set>
#include "gridBoundary.hh"
#include "../Math/vectorMath.hh"

#define IMPLEMENT_APPLY_INTERFACE(type, TypeName, T) \
    void Apply##TypeName(type* field) override { ApplyThis<T>(field); }

// Base class for Grid Boundary
template <int dim>
class ReflectingGridBoundary : public GridBoundary<dim> {
public:
    using Vector      = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;
    using Complex     = std::complex<double>;
    using ComplexField= Field<Complex>;

private:
    // Face-based reflecting BC (domain boundaries)
    std::vector<std::vector<int>> boundaryLists;
    std::vector<std::vector<int>> interiorLists;

    // Interior Neumann obstacle support
    std::vector<int> obstacleIds;
    std::vector<std::vector<int>> obstacleNeighbors;  // exterior neighbors per obstacle cell
    bool initialized = false;

    template <typename T>
    void ApplyThis(Field<T>* field) {
        // ── face-based reflecting BC ───────────────────────────────────────
        #pragma omp parallel for
        for (int i = 0; i < (int)boundaryLists.size(); ++i) {
            if (!this->activeFaces[i]) continue;
            const auto& b    = boundaryLists[i];
            const auto& iids = interiorLists[i];
            for (size_t j = 0; j < b.size(); ++j) {
                if constexpr (std::is_same_v<T, Vector>) {
                    Vector val = field->getValue(iids[j]);
                    val[i/2] *= -1;  // flip normal component
                    field->setValue(b[j], val);
                } else if constexpr (std::is_same_v<T, double> ||
                                     std::is_same_v<T, Complex>) {
                    field->setValue(b[j], field->getValue(iids[j]));
                }
            }
        }

        // ── interior obstacle BC ────────────────────────────────────────────
        // Scalars (e.g. phi/xi for the wave equation, or density/specificInternalEnergy
        // for hydro): each obstacle cell is set to the average of its exterior
        // (non-obstacle) neighbors, giving ∂φ/∂n = 0 at the obstacle surface —
        // the rigid-wall condition.
        //
        // Vectors (e.g. velocity for hydro): held at zero instead of averaged.
        // Averaging would just reproduce whatever the surrounding flow is doing
        // (in a uniform freestream, "average of neighbors" == freestream), which
        // makes the obstacle invisible to the flow. Pinning velocity to zero
        // gives a stationary rigid body that the surrounding fluid must actually
        // flow around.
        if (!initialized) return;
        #pragma omp parallel for
        for (int i = 0; i < (int)obstacleIds.size(); ++i) {
            if constexpr (std::is_same_v<T, Vector>) {
                field->setValue(obstacleIds[i], Vector());
            } else {
                const auto& nbrs = obstacleNeighbors[i];
                if (nbrs.empty()) continue;
                T sum{};
                for (int n : nbrs) sum = sum + field->getValue(n);
                field->setValue(obstacleIds[i], sum * (1.0 / nbrs.size()));
            }
        }
    }

public:
    ReflectingGridBoundary(Mesh::Grid<dim>* grid) :
        GridBoundary<dim>(grid) {
        if constexpr (dim == 1) {
            boundaryLists.push_back(grid->leftMost());
            interiorLists.push_back({1});

            boundaryLists.push_back(grid->rightMost());
            interiorLists.push_back({grid->rightMost()[0]-1});
        }
        else if constexpr (dim == 2) {
            std::vector<int> left   = grid->leftMost();
            std::vector<int> right  = grid->rightMost();
            std::vector<int> top    = grid->topMost();
            std::vector<int> bottom = grid->bottomMost();

            std::vector<int> leftInt, rightInt, topInt, bottomInt;

            for (int j = 0; j < (int)right.size(); ++j) {
                leftInt.push_back(grid->index(1, j));
                rightInt.push_back(grid->index(grid->size_x() - 2, j));
            }

            for (int i = 0; i < (int)bottom.size(); ++i) {
                topInt.push_back(grid->index(i, 1));
                bottomInt.push_back(grid->index(i, grid->size_y() - 2));
            }

            boundaryLists = {left, right, top, bottom};
            interiorLists = {leftInt, rightInt, topInt, bottomInt};
        }
        else if constexpr (dim == 3) {
            std::vector<int> left   = grid->leftMost();
            std::vector<int> right  = grid->rightMost();
            std::vector<int> top    = grid->topMost();
            std::vector<int> bottom = grid->bottomMost();
            std::vector<int> front  = grid->frontMost();
            std::vector<int> back   = grid->backMost();

            std::vector<int> leftInt, rightInt, topInt, bottomInt, frontInt, backInt;

            int Nx = grid->size_x();
            int Ny = grid->size_y();
            int Nz = grid->size_z();

            // Iteration order here must match Grid<3>::findBoundaries exactly,
            // since ApplyThis pairs boundaryLists[i][j] with interiorLists[i][j]
            // by position: j outer / k inner for left-right (matches lm/rm).
            for (int j = 0; j < Ny; ++j) {
                for (int k = 0; k < Nz; ++k) {
                    leftInt.push_back(grid->index(1, j, k));
                    rightInt.push_back(grid->index(Nx - 2, j, k));
                }
            }

            // i outer / k inner for top-bottom (matches tm/bm).
            for (int i = 0; i < Nx; ++i) {
                for (int k = 0; k < Nz; ++k) {
                    topInt.push_back(grid->index(i, 1, k));
                    bottomInt.push_back(grid->index(i, Ny - 2, k));
                }
            }

            // i outer / j inner for front-back (matches fm/km).
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    frontInt.push_back(grid->index(i, j, 1));
                    backInt.push_back(grid->index(i, j, Nz - 2));
                }
            }

            boundaryLists = {left, right, top, bottom, front, back};
            interiorLists = {leftInt, rightInt, topInt, bottomInt, frontInt, backInt};
        }
    }

    virtual ~ReflectingGridBoundary() = default;

    // Returns the obstacle cell IDs so the physics update loop can exclude them.
    virtual std::vector<int> getObstacleIds() const override { return obstacleIds; }

    // Called once before the first step.  Builds the exterior-neighbor lookup
    // table for all registered obstacle cells.
    virtual void ZeroTimeInitialize(NodeList* nodeList) override {
        if (obstacleIds.empty()) return;

        // Deduplicate obstacle IDs (multiple addBox calls may overlap).
        std::unordered_set<int> obstacleSet(obstacleIds.begin(), obstacleIds.end());
        obstacleIds.assign(obstacleSet.begin(), obstacleSet.end());

        // For each obstacle cell, collect neighbors that are NOT obstacles.
        obstacleNeighbors.resize(obstacleIds.size());
        for (int i = 0; i < (int)obstacleIds.size(); ++i) {
            obstacleNeighbors[i].clear();
            for (int n : this->grid->getNeighboringCells(obstacleIds[i])) {
                if (n < 0) continue;  // off the domain edge, not a real neighbor
                if (obstacleSet.count(n) == 0)
                    obstacleNeighbors[i].push_back(n);
            }
        }
        initialized = true;
    }

    // Register an axis-aligned box of cells as a Neumann (reflecting) obstacle.
    // p1/p2 are in real-space coordinates (cell centres at i*dx+0.5*dx, etc.).
    virtual void addBox(Vector p1, Vector p2) {
        for (int idx = 0; idx < this->grid->size(); ++idx) {
            Vector pos = this->grid->getPosition(idx);
            bool inside = true;
            for (int d = 0; d < dim; ++d)
                if (pos[d] < p1[d] || pos[d] > p2[d]) { inside = false; break; }
            if (inside) obstacleIds.push_back(idx);
        }
    }

    virtual void removeBox(Vector p1, Vector p2) {
        std::vector<int> toRemove;
        for (int idx = 0; idx < this->grid->size(); ++idx) {
            Vector pos = this->grid->getPosition(idx);
            bool inside = true;
            for (int d = 0; d < dim; ++d)
                if (pos[d] < p1[d] || pos[d] > p2[d]) { inside = false; break; }
            if (inside) toRemove.push_back(idx);
        }
        for (int r : toRemove) {
            auto it = std::find(obstacleIds.begin(), obstacleIds.end(), r);
            if (it != obstacleIds.end()) obstacleIds.erase(it);
        }
    }

    virtual void addSphere(Vector p, double radius) {
        for (int idx = 0; idx < this->grid->size(); ++idx) {
            Vector pos = this->grid->getPosition(idx);
            if ((pos - p).mag2() <= radius * radius)
                obstacleIds.push_back(idx);
        }
    }

    virtual void removeSphere(Vector p, double radius) {
        std::vector<int> toRemove;
        for (int idx = 0; idx < this->grid->size(); ++idx) {
            Vector pos = this->grid->getPosition(idx);
            if ((pos - p).mag2() <= radius * radius)
                toRemove.push_back(idx);
        }
        for (int r : toRemove) {
            auto it = std::find(obstacleIds.begin(), obstacleIds.end(), r);
            if (it != obstacleIds.end()) obstacleIds.erase(it);
        }
    }

    IMPLEMENT_APPLY_INTERFACE(ScalarField, Scalar, double)
    IMPLEMENT_APPLY_INTERFACE(VectorField, Vector, Vector)
    IMPLEMENT_APPLY_INTERFACE(ComplexField, Complex, Complex)
};
