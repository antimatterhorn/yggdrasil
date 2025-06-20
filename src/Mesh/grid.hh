// Copyright (C) 2025  Cody Raskin

#pragma once

#include <vector>
#include <string>
#include "../Math/vectorMath.hh"
#include "../DataBase/field.hh"
#include "../DataBase/nodeList.hh"

namespace Mesh {
    template <int dim>
    class Grid {
    private:
        std::vector<std::shared_ptr<FieldBase>> _extraFields;
        std::vector<int> lm,rm,tm,bm,fm,km;
    public:
        using Vector = Lin::Vector<dim>;
        int nx; // Number of grid cells in x-direction
        int ny; // Number of grid cells in y-direction
        int nz; // Number of grid cells in z-direction

        double dx; // Grid spacing in x-direction
        double dy; // Grid spacing in y-direction
        double dz; // Grid spacing in z-direction

        Field<Vector> gridPositions;

        Grid(int num_cells_x, double spacing_x);
        Grid(int num_cells_x, int num_cells_y, double spacing_x, double spacing_y);
        Grid(int num_cells_x, int num_cells_y, int num_cells_z, double spacing_x, double spacing_y, double spacing_z);

        void initializeGrid();
        void setOrigin(Lin::Vector<dim> origin);

        int index(int i, int j = 0, int k = 0) const;

        inline int getnx() const    { return nx; };
        inline int getny() const    { return ny; };
        inline int getnz() const    { return nz; };
        inline int size_x() const   { return nx; };
        inline int size_y() const   { return ny; };
        inline int size_z() const   { return nz; };
        inline int size() const     { return nx*ny*nz; };
        inline double getdx() const { return dx; };
        inline double getdy() const { return dy; };
        inline double getdz() const { return dz; };
        inline std::vector<int> leftMost() const    { return lm; };
        inline std::vector<int> rightMost() const   { return rm; };
        inline std::vector<int> topMost() const     { return tm; };
        inline std::vector<int> bottomMost() const  { return bm; };
        inline std::vector<int> frontMost() const   { return fm; };
        inline std::vector<int> backMost() const    { return km; };

        inline double cellVolume(int idx) const { return dx * dy * dz; };  // for now, grids are homogeneous
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

        double spacing(int axis) const;
        Lin::Vector<dim> getPosition(int id);

        std::vector<int> getNeighboringCells(int idx) const;
        inline std::vector<int> neighbors(int idx ) const { return getNeighboringCells(idx); };
        std::array<int, 3> indexToCoordinates(int idx) const;

        void findBoundaries(const int buffer);
        bool onBoundary(const int idx);
        void assignPositions(NodeList* nodeList);

        template <typename T>
        void insertField(const std::string& name);

        template <typename T>
        Field<T>* getField(const std::string& name);

        template <typename T>
        T laplacian(int idx, Field<T>* field) const;

        Vector gradient(int idx, Field<double>* field) const {
            Vector grad;

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
    };
}


#include "grid.cc"
