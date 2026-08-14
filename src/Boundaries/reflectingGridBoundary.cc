// Copyright (C) 2026  Cody Raskin

#pragma once
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
        // boundaryLists[i] are the ghost cells for face i; interiorLists[i][j]
        // is boundaryLists[i][j]'s mirror source -- the real first/last
        // logical cell along that axis (index 0 / n-1).
        if constexpr (dim == 1) {
            boundaryLists.push_back(grid->lowMost(0));
            interiorLists.push_back({grid->index(0)});

            boundaryLists.push_back(grid->highMost(0));
            interiorLists.push_back({grid->index(grid->size_x() - 1)});
        }
        else if constexpr (dim == 2) {
            std::vector<int> xLow  = grid->lowMost(0);
            std::vector<int> xHigh = grid->highMost(0);
            std::vector<int> yLow  = grid->lowMost(1);
            std::vector<int> yHigh = grid->highMost(1);

            std::vector<int> xLowInt, xHighInt, yLowInt, yHighInt;

            int Nx = grid->size_x();
            int Ny = grid->size_y();

            for (int j = 0; j < Ny; ++j) {
                xLowInt.push_back(grid->index(0, j));
                xHighInt.push_back(grid->index(Nx - 1, j));
            }

            for (int i = 0; i < Nx; ++i) {
                yLowInt.push_back(grid->index(i, 0));
                yHighInt.push_back(grid->index(i, Ny - 1));
            }

            boundaryLists = {xLow, xHigh, yLow, yHigh};
            interiorLists = {xLowInt, xHighInt, yLowInt, yHighInt};
        }
        else if constexpr (dim == 3) {
            std::vector<int> xLow  = grid->lowMost(0);
            std::vector<int> xHigh = grid->highMost(0);
            std::vector<int> yLow  = grid->lowMost(1);
            std::vector<int> yHigh = grid->highMost(1);
            std::vector<int> zLow  = grid->lowMost(2);
            std::vector<int> zHigh = grid->highMost(2);

            std::vector<int> xLowInt, xHighInt, yLowInt, yHighInt, zLowInt, zHighInt;

            int Nx = grid->size_x();
            int Ny = grid->size_y();
            int Nz = grid->size_z();

            // Iteration order here must match Grid<3>::buildGhostLists exactly,
            // since ApplyThis pairs boundaryLists[i][j] with interiorLists[i][j]
            // by position: j outer / k inner for the x faces.
            for (int j = 0; j < Ny; ++j) {
                for (int k = 0; k < Nz; ++k) {
                    xLowInt.push_back(grid->index(0, j, k));
                    xHighInt.push_back(grid->index(Nx - 1, j, k));
                }
            }

            // i outer / k inner for the y faces.
            for (int i = 0; i < Nx; ++i) {
                for (int k = 0; k < Nz; ++k) {
                    yLowInt.push_back(grid->index(i, 0, k));
                    yHighInt.push_back(grid->index(i, Ny - 1, k));
                }
            }

            // i outer / j inner for the z faces.
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    zLowInt.push_back(grid->index(i, j, 0));
                    zHighInt.push_back(grid->index(i, j, Nz - 1));
                }
            }

            boundaryLists = {xLow, xHigh, yLow, yHigh, zLow, zHigh};
            interiorLists = {xLowInt, xHighInt, yLowInt, yHighInt, zLowInt, zHighInt};
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
    // Restricted to logical cells -- an obstacle is carved out of real domain
    // interior, never out of the ghost halo (that's what setFaces is for).
    virtual void addBox(Vector p1, Vector p2) {
        for (int idx = 0; idx < this->grid->size(); ++idx) {
            if (this->grid->isGhost(idx)) continue;
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
            if (this->grid->isGhost(idx)) continue;
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
            if (this->grid->isGhost(idx)) continue;
            Vector pos = this->grid->getPosition(idx);
            if ((pos - p).mag2() <= radius * radius)
                obstacleIds.push_back(idx);
        }
    }

    virtual void removeSphere(Vector p, double radius) {
        std::vector<int> toRemove;
        for (int idx = 0; idx < this->grid->size(); ++idx) {
            if (this->grid->isGhost(idx)) continue;
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
