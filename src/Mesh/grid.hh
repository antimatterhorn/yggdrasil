// Copyright (C) 2026  Cody Raskin

#pragma once

#include <vector>
#include <string>
#include "../Math/vectorMath.hh"
#include "../DataBase/field.hh"
#include "../DataBase/nodeList.hh"
#include "geometry.hh"

namespace Mesh {

    template <int dim>
    class Grid {
    private:
        std::vector<std::shared_ptr<FieldBase>> _extraFields;
        std::array<std::vector<int>, 3> lowIds, highIds;
        Geometry _geometry = Geometry::Cartesian;
        // Allocated (logical + ghost) extents. nx/ny/nz below stay the
        // user-requested logical extent; these are nx+2*ghostWidth etc. on
        // every axis that exists for this dim, unpadded (== nx/ny/nz) otherwise.
        int nxTotal = 0, nyTotal = 0, nzTotal = 0;
    public:
        using Vector = Lin::Vector<dim>;
        int nx; // Number of logical (real, user-visible) grid cells in x
        int ny; // Number of logical grid cells in y
        int nz; // Number of logical grid cells in z

        double dx; // Grid spacing in x-direction
        double dy; // Grid spacing in y-direction
        double dz; // Grid spacing in z-direction

        // Depth of the ghost halo surrounding the logical domain on every
        // side. A domain-edge face has only one real neighbor beyond it, so
        // GridHydroHLLC/KT's 4-cell stencil falls back to piecewise-constant
        // at that face (see their computeFlux).
        static constexpr int ghostWidth = 1;

        Field<Vector> gridPositions;

        Grid(int num_cells_x, double spacing_x);
        Grid(int num_cells_x, int num_cells_y, double spacing_x, double spacing_y,
             Geometry geometry = Geometry::Cartesian);
        Grid(int num_cells_x, int num_cells_y, int num_cells_z, double spacing_x, double spacing_y, double spacing_z);

        void initializeGrid();
        void setOrigin(Lin::Vector<dim> origin);

        inline Geometry geometry() const { return _geometry; }

        // Logical coordinates: i in [-ghostWidth, nx-1+ghostWidth] (and likewise
        // j, k on any axis this dim has), 0..n-1 being the real domain and the
        // rest being ghost cells. Maps to the flat storage index.
        int index(int i, int j = 0, int k = 0) const;

        inline int getnx() const    { return nx; };
        inline int getny() const    { return ny; };
        inline int getnz() const    { return nz; };
        inline int size_x() const   { return nx; };
        inline int size_y() const   { return ny; };
        inline int size_z() const   { return nz; };
        // Total allocated storage (logical cells + ghost halo) -- what a
        // NodeList living on this grid should be sized to, e.g.
        // NodeList(grid.size()).
        inline int size() const     { return nxTotal * nyTotal * nzTotal; };
        // Real (non-ghost) cell count -- nx*ny*nz.
        inline int logicalSize() const { return nx * ny * nz; };
        inline double getdx() const { return dx; };
        inline double getdy() const { return dy; };
        inline double getdz() const { return dz; };
        // Ghost cell indices one halo layer out from the low/high side of
        // `axis` (0=x, 1=y, 2=z), matching getNeighboringCells' [low-x, high-x,
        // low-y, ...] ordering.
        inline std::vector<int> lowMost(int axis) const  { return lowIds[axis]; };
        inline std::vector<int> highMost(int axis) const { return highIds[axis]; };

        // Cell-centre radius (the axis-0 coordinate), used as r by the RZ metric.
        inline double cellRadius(int idx) const { return gridPositions[idx][0]; }

        // Per-radian annular volume r*dr*dz in RZ (the 2*pi cancels in the
        // flux/V update, so it is omitted); uniform dx*dy*dz otherwise.
        inline double cellVolume(int idx) const {
            if constexpr (dim == 2) {
                if (_geometry == Geometry::CylindricalRZ)
                    return cellRadius(idx) * dx * dy;
            }
            return dx * dy * dz;
        }

        // Area of the face shared by adjacent cells i and j (orientation inferred
        // from their index offset).
        inline double faceArea(int i, int j) const {
            auto ci = indexToCoordinates(i);
            auto cj = indexToCoordinates(j);

            int dx_ = std::abs(ci[0] - cj[0]);
            int dy_ = std::abs(ci[1] - cj[1]);
            int dz_ = std::abs(ci[2] - cj[2]);

            if constexpr (dim == 1) {
                return 1.0;  // face area is unitless in 1D
            }
            else if constexpr (dim == 2) {
                if (_geometry == Geometry::CylindricalRZ) {
                    if (dx_ == 1 && dy_ == 0)   // r-face between the two cells
                        return 0.5 * (cellRadius(i) + cellRadius(j)) * dy;
                    if (dy_ == 1 && dx_ == 0)   // z-face (annulus of the shared cell)
                        return cellRadius(i) * dx;
                }
                if (dx_ == 1 && dy_ == 0) return dy;
                if (dy_ == 1 && dx_ == 0) return dx;
            }
            else if constexpr (dim == 3) {
                if (dx_ == 1 && dy_ == 0 && dz_ == 0) return dy * dz;
                if (dx_ == 0 && dy_ == 1 && dz_ == 0) return dx * dz;
                if (dx_ == 0 && dy_ == 0 && dz_ == 1) return dx * dy;
            }

            // Diagonal or non-adjacent — undefined
            return 0.0;
        }

        // Area of the face between cells iL and iR normal to `axis`.
        inline double faceArea(int iL, int iR, int axis) const {
            if constexpr (dim == 1) {
                return 1.0;
            }
            else if constexpr (dim == 2) {
                if (_geometry == Geometry::CylindricalRZ) {
                    if (axis == 0) return 0.5 * (cellRadius(iL) + cellRadius(iR)) * dy;  // r-face
                    if (axis == 1) return cellRadius(iL) * dx;                           // z-face
                }
                return (axis == 0) ? dy : dx;
            }
            else if constexpr (dim == 3) {
                if (axis == 0) return dy * dz;
                if (axis == 1) return dx * dz;
                return dx * dy;
            }
            return 0.0;
        }

        double spacing(int axis) const;
        Lin::Vector<dim> getPosition(int id);

        std::array<int, 2*dim> getNeighboringCells(int idx) const;
        inline std::array<int, 2*dim> neighbors(int idx ) const { return getNeighboringCells(idx); };
        std::array<int, 3> indexToCoordinates(int idx) const;

        void buildGhostLists();
        // Whether cell `idx` lies outside the logical [0,nx)x[0,ny)x[0,nz)
        // domain -- i.e. is part of the ghost halo rather than real state.
        bool isGhost(int idx) const;
        inline bool isInterior(int idx) const { return !isGhost(idx); }
        void assignPositions(NodeList* nodeList);

        template <typename T>
        void insertField(const std::string& name);

        template <typename T>
        Field<T>* getField(const std::string& name);

        template <typename T>
        T laplacian(int idx, Field<T>* field) const;

        template <typename T, typename F>
        T laplacian(int idx, F&& valueAt) const;

        Vector gradient(int idx, Field<double>* field) const;

        template <typename F>
        Vector gradient(int idx, F&& valueAt) const;
    };
}


#include "grid.cc"
