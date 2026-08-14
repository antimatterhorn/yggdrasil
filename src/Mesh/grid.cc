// Copyright (C) 2026  Cody Raskin

#ifndef GRID_CC
#define GRID_CC

#include "grid.hh"

namespace Mesh {
    template <int dim>
    Grid<dim>::Grid(int num_cells_x, double spacing_x)
        : nx(num_cells_x), ny(1), nz(1), dx(spacing_x), dy(1), dz(1) {
        initializeGrid();
    }

    template <int dim>
    Grid<dim>::Grid(int num_cells_x, int num_cells_y, double spacing_x, double spacing_y,
                    Geometry geometry)
        : _geometry(geometry), nx(num_cells_x), ny(num_cells_y), nz(1), dx(spacing_x), dy(spacing_y), dz(1) {
        initializeGrid();
    }

    template <int dim>
    Grid<dim>::Grid(int num_cells_x, int num_cells_y, int num_cells_z, double spacing_x, double spacing_y, double spacing_z)
        : nx(num_cells_x), ny(num_cells_y), nz(num_cells_z), dx(spacing_x), dy(spacing_y), dz(spacing_z) {
        initializeGrid();
    }

    template <int dim>
    void
    Grid<dim>::initializeGrid() {
        nxTotal = nx + 2 * ghostWidth;
        nyTotal = (dim >= 2) ? ny + 2 * ghostWidth : ny;
        nzTotal = (dim >= 3) ? nz + 2 * ghostWidth : nz;

        gridPositions = Field<Lin::Vector<dim>>("gridPosition");

        // Cell-centre positions over the full padded (logical + ghost) domain,
        // in the same k-major/j/i-minor flat order as index(), so ghost cells
        // get sensible extrapolated positions (the i*dx+0.5*dx formula works
        // unchanged for negative i) without any special-casing.
        const int kLo = (dim >= 3) ? -ghostWidth : 0, kHi = (dim >= 3) ? nz - 1 + ghostWidth : 0;
        const int jLo = (dim >= 2) ? -ghostWidth : 0, jHi = (dim >= 2) ? ny - 1 + ghostWidth : 0;
        const int iLo = -ghostWidth, iHi = nx - 1 + ghostWidth;

        for (int k = kLo; k <= kHi; ++k) {
            for (int j = jLo; j <= jHi; ++j) {
                for (int i = iLo; i <= iHi; ++i) {
                    Lin::Vector<dim> gridPosition;
                    for (int d = 0; d < dim; ++d) {
                        gridPosition.values[d] = (d == 0 ? i * dx + 0.5 * dx : (d == 1 ? j * dy + 0.5 * dy : k * dz + 0.5 * dz));
                    }
                    gridPositions.addValue(gridPosition);
                }
            }
        }

        buildGhostLists();
    }

    template <int dim>
    void 
    Grid<dim>::setOrigin(Lin::Vector<dim> origin) {
        #pragma omp parallel for
        for (int i = 0; i < gridPositions.getSize(); ++i) {
            gridPositions[i] -= origin;
        }
    }

    // Specialization for 1D
    template <>
    int Grid<1>::index(int i, int, int) const {
        return i + ghostWidth;
    }

    // Specialization for 2D
    template <>
    int Grid<2>::index(int i, int j, int) const {
        return (j + ghostWidth) * nxTotal + (i + ghostWidth);
    }

    // Specialization for 3D
    template <>
    int Grid<3>::index(int i, int j, int k) const {
        return ((k + ghostWidth) * nyTotal + (j + ghostWidth)) * nxTotal + (i + ghostWidth);
    }

    // Generic template definition (if you need for other dimensions)
    template <int dim>
    int Grid<dim>::index(int i, int j, int k) const {
        // This is a generic case. You can define custom behavior or throw an error if needed.
        static_assert(dim <= 3, "Only up to 3 dimensions are supported");
        return -1; // Error value
    }

    template <int dim>
    double 
    Grid<dim>::spacing(int axis) const {
        if (axis==0)
            return dx;
        else if (axis==1)
            return dy;
        else if (axis==2)
            return dz;
        return 1;
    }

    template <int dim>
    Lin::Vector<dim> 
    Grid<dim>::getPosition(int id) {
        return gridPositions[id];
    }

    // Always returns exactly 2*dim entries, in fixed axis-major order
    // [-axis0, +axis0, -axis1, +axis1, -axis2, +axis2], using -1 for any
    // direction that runs off the padded (logical + ghost) array. Callers
    // (e.g. GridHydroKT's slopeLimitedValue) rely on both the fixed
    // size/position and the -1 sentinel to detect a missing neighbor at the
    // edge of the ghost halo -- do not go back to conditionally omitting
    // entries, which silently reindexes every neighbor after the gap and
    // invites out-of-bounds reads.
    template <int dim>
    std::array<int, 2*dim>
    Grid<dim>::getNeighboringCells(int idx) const {
        std::array<int, 3> coords = indexToCoordinates(idx);
        std::array<int, 2*dim> neighbors;
        neighbors.fill(-1);

        neighbors[0] = (coords[0] - 1 >= -ghostWidth)      ? index(coords[0] - 1, coords[1], coords[2]) : -1;
        neighbors[1] = (coords[0] + 1 <= nx - 1 + ghostWidth) ? index(coords[0] + 1, coords[1], coords[2]) : -1;

        if constexpr (dim > 1) {
            neighbors[2] = (coords[1] - 1 >= -ghostWidth)      ? index(coords[0], coords[1] - 1, coords[2]) : -1;
            neighbors[3] = (coords[1] + 1 <= ny - 1 + ghostWidth) ? index(coords[0], coords[1] + 1, coords[2]) : -1;
        }

        if constexpr (dim == 3) {
            neighbors[4] = (coords[2] - 1 >= -ghostWidth)      ? index(coords[0], coords[1], coords[2] - 1) : -1;
            neighbors[5] = (coords[2] + 1 <= nz - 1 + ghostWidth) ? index(coords[0], coords[1], coords[2] + 1) : -1;
        }

        return neighbors;
    }

    // Inverse of index(): flat storage index -> logical (i,j,k). Ghost cells
    // decode to coordinates outside [0,n) on the corresponding axis.
    template <int dim>
    std::array<int, 3>
    Grid<dim>::indexToCoordinates(int idx) const {
        std::array<int, 3> coords;
        coords.fill(0);
        if constexpr (dim == 3) {
            coords[2] = idx / (nxTotal * nyTotal); // Compute k
            idx -= coords[2] * nxTotal * nyTotal;
            coords[2] -= ghostWidth;
        }
        if constexpr (dim > 1)
            coords[1] = idx / nxTotal - ghostWidth; // Compute j
        coords[0] = idx % nxTotal - ghostWidth; // Compute i
        return coords;
    }

    // Rebuilds the per-axis ghost-cell index lists (lowMost/highMost) from the
    // current logical extents. One layer per call, `ghostWidth` layers deep.
    template <int dim>
    void
    Grid<dim>::buildGhostLists() {
        for (int d = 0; d < 3; ++d) {
            lowIds[d].clear();
            highIds[d].clear();
        }

        if constexpr (dim == 1) {
            for (int b = 0; b < ghostWidth; ++b) {
                lowIds[0].push_back(index(-1 - b));
                highIds[0].push_back(index(nx + b));
            }
        }
        else if constexpr (dim == 2) {
            for (int b = 0; b < ghostWidth; ++b) {
                for (int j = 0; j < ny; ++j) {
                    lowIds[0].push_back(index(-1 - b, j));
                    highIds[0].push_back(index(nx + b, j));
                }
                for (int i = 0; i < nx; ++i) {
                    lowIds[1].push_back(index(i, -1 - b));
                    highIds[1].push_back(index(i, ny + b));
                }
            }
        }
        else if constexpr (dim == 3) {
            for (int b = 0; b < ghostWidth; ++b) {
                for (int j = 0; j < ny; ++j) {
                    for (int k = 0; k < nz; ++k) {
                        lowIds[0].push_back(index(-1 - b, j, k));
                        highIds[0].push_back(index(nx + b, j, k));
                    }
                }
                for (int i = 0; i < nx; ++i) {
                    for (int k = 0; k < nz; ++k) {
                        lowIds[1].push_back(index(i, -1 - b, k));
                        highIds[1].push_back(index(i, ny + b, k));
                    }
                }
                for (int i = 0; i < nx; ++i) {
                    for (int j = 0; j < ny; ++j) {
                        lowIds[2].push_back(index(i, j, -1 - b));
                        highIds[2].push_back(index(i, j, nz + b));
                    }
                }
            }
        }
    }

    template <int dim>
    bool
    Grid<dim>::isGhost(int idx) const {
        auto c = indexToCoordinates(idx);
        if (c[0] < 0 || c[0] >= nx) return true;
        if constexpr (dim > 1) if (c[1] < 0 || c[1] >= ny) return true;
        if constexpr (dim > 2) if (c[2] < 0 || c[2] >= nz) return true;
        return false;
    }

    template <int dim>
    void
    Grid<dim>::assignPositions(NodeList* nodeList) {
        using Vector = Lin::Vector<dim>;
        using VectorField = Field<Vector>;
        if (nodeList->getField<Vector>("position") == nullptr)
            nodeList->insertField<Vector>("position");
                
        VectorField* position = nodeList->getField<Vector>("position");
        #pragma omp parallel for
        for (int i = 0; i < position->size(); ++i) {
            position->setValue(i, this->getPosition(i));
        }
    }

    template <int dim>
    template <typename T>
    void Grid<dim>::insertField(const std::string& name) {
        auto newField = std::make_shared<Field<T>>(name, this->size());
        _extraFields.push_back(newField);
    }

    template <int dim>
    template <typename T>
    Field<T>* Grid<dim>::getField(const std::string& name) {
        for (const auto& field : _extraFields) {
            auto specificField = std::dynamic_pointer_cast<Field<T>>(field);
            if (specificField && specificField->getName() == name) {
                return specificField.get();
            }
        }
        return nullptr;
    }

    template <int dim>
    template <typename T>
    T Grid<dim>::laplacian(int idx, Field<T>* field) const {
        if constexpr (dim == 1) {
            double dx = this->dx;
            int left  = idx - 1;
            int right = idx + 1;
            return (field->getValue(right) - 2.0 * field->getValue(idx) + field->getValue(left)) / (dx * dx);
        }

        if constexpr (dim == 2) {
            auto [ix, iy, _] = this->indexToCoordinates(idx);
            double dx = this->dx;
            double dy = this->dy;

            int iL = this->index(ix - 1, iy);
            int iR = this->index(ix + 1, iy);
            int iB = this->index(ix, iy - 1);
            int iT = this->index(ix, iy + 1);

            T center = field->getValue(idx);
            return (field->getValue(iR) - 2.0 * center + field->getValue(iL)) / (dx * dx) +
                (field->getValue(iT) - 2.0 * center + field->getValue(iB)) / (dy * dy);
        }

        if constexpr (dim == 3) {
            auto [ix, iy, iz] = this->indexToCoordinates(idx);
            double dx = this->dx;
            double dy = this->dy;
            double dz = this->dz;

            int iL = this->index(ix - 1, iy, iz);
            int iR = this->index(ix + 1, iy, iz);
            int iB = this->index(ix, iy - 1, iz);
            int iT = this->index(ix, iy + 1, iz);
            int iD = this->index(ix, iy, iz - 1);
            int iU = this->index(ix, iy, iz + 1);

            T center = field->getValue(idx);
            return (field->getValue(iR) - 2.0 * center + field->getValue(iL)) / (dx * dx) +
                (field->getValue(iT) - 2.0 * center + field->getValue(iB)) / (dy * dy) +
                (field->getValue(iU) - 2.0 * center + field->getValue(iD)) / (dz * dz);
        }
    }

    template <int dim>
    template <typename T, typename F>
    T Grid<dim>::laplacian(int idx, F&& valueAt) const {
        if constexpr (dim == 1) {
            double dx = this->dx;
            int left  = idx - 1;
            int right = idx + 1;
            return (valueAt(right) - 2.0 * valueAt(idx) + valueAt(left)) / (dx * dx);
        }

        if constexpr (dim == 2) {
            auto [ix, iy, _] = this->indexToCoordinates(idx);
            double dx = this->dx;
            double dy = this->dy;

            int iL = this->index(ix - 1, iy);
            int iR = this->index(ix + 1, iy);
            int iB = this->index(ix, iy - 1);
            int iT = this->index(ix, iy + 1);

            return (valueAt(iR) - 2.0 * valueAt(idx) + valueAt(iL)) / (dx * dx) +
                (valueAt(iT) - 2.0 * valueAt(idx) + valueAt(iB)) / (dy * dy);
        }

        if constexpr (dim == 3) {
            auto [ix, iy, iz] = this->indexToCoordinates(idx);
            double dx = this->dx;
            double dy = this->dy;
            double dz = this->dz;

            int iL = this->index(ix - 1, iy, iz);
            int iR = this->index(ix + 1, iy, iz);
            int iB = this->index(ix, iy - 1, iz);
            int iT = this->index(ix, iy + 1, iz);
            int iD = this->index(ix, iy, iz - 1);
            int iU = this->index(ix, iy, iz + 1);

            return (valueAt(iR) - 2.0 * valueAt(idx) + valueAt(iL)) / (dx * dx) +
                (valueAt(iT) - 2.0 * valueAt(idx) + valueAt(iB)) / (dy * dy) +
                (valueAt(iU) - 2.0 * valueAt(idx) + valueAt(iD)) / (dz * dz);
        }
    }

    template <int dim> Lin::Vector<dim> 
    Grid<dim>::gradient(int idx, Field<double>* field) const {
        Lin::Vector<dim> grad;

        if constexpr (dim == 1) {
            double dx = this->dx;
            int left = idx - 1;
            int right = idx + 1;
            grad[0] = (field->getValue(right) - field->getValue(left)) / (2.0 * dx);
        }

        if constexpr (dim == 2) {
            auto [ix, iy, _] = this->indexToCoordinates(idx);
            double dx = this->dx;
            double dy = this->dy;

            int iL = this->index(ix - 1, iy);
            int iR = this->index(ix + 1, iy);
            int iB = this->index(ix, iy - 1);
            int iT = this->index(ix, iy + 1);

            grad[0] = (field->getValue(iR) - field->getValue(iL)) / (2.0 * dx);
            grad[1] = (field->getValue(iT) - field->getValue(iB)) / (2.0 * dy);
        }

        if constexpr (dim == 3) {
            auto [ix, iy, iz] = this->indexToCoordinates(idx);
            double dx = this->dx;
            double dy = this->dy;
            double dz = this->dz;

            int iL = this->index(ix - 1, iy, iz);
            int iR = this->index(ix + 1, iy, iz);
            int iB = this->index(ix, iy - 1, iz);
            int iT = this->index(ix, iy + 1, iz);
            int iD = this->index(ix, iy, iz - 1);
            int iU = this->index(ix, iy, iz + 1);

            grad[0] = (field->getValue(iR) - field->getValue(iL)) / (2.0 * dx);
            grad[1] = (field->getValue(iT) - field->getValue(iB)) / (2.0 * dy);
            grad[2] = (field->getValue(iU) - field->getValue(iD)) / (2.0 * dz);
        }

        return grad;
    }

    template <int dim>
    template <typename F>
    Lin::Vector<dim> Grid<dim>::gradient(int idx, F&& valueAt) const {
        Lin::Vector<dim> grad;

        if constexpr (dim == 1) {
            double dx = this->dx;
            int left  = idx - 1;
            int right = idx + 1;
            grad[0] = (valueAt(right) - valueAt(left)) / (2.0 * dx);
        }

        if constexpr (dim == 2) {
            auto [ix, iy, _] = this->indexToCoordinates(idx);
            double dx = this->dx;
            double dy = this->dy;

            int iL = this->index(ix - 1, iy);
            int iR = this->index(ix + 1, iy);
            int iB = this->index(ix, iy - 1);
            int iT = this->index(ix, iy + 1);

            grad[0] = (valueAt(iR) - valueAt(iL)) / (2.0 * dx);
            grad[1] = (valueAt(iT) - valueAt(iB)) / (2.0 * dy);
        }

        if constexpr (dim == 3) {
            auto [ix, iy, iz] = this->indexToCoordinates(idx);
            double dx = this->dx;
            double dy = this->dy;
            double dz = this->dz;

            int iL = this->index(ix - 1, iy, iz);
            int iR = this->index(ix + 1, iy, iz);
            int iB = this->index(ix, iy - 1, iz);
            int iT = this->index(ix, iy + 1, iz);
            int iD = this->index(ix, iy, iz - 1);
            int iU = this->index(ix, iy, iz + 1);

            grad[0] = (valueAt(iR) - valueAt(iL)) / (2.0 * dx);
            grad[1] = (valueAt(iT) - valueAt(iB)) / (2.0 * dy);
            grad[2] = (valueAt(iU) - valueAt(iD)) / (2.0 * dz);
        }

        return grad;
    }


}


#endif