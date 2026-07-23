// Copyright (C) 2026  Cody Raskin

#ifndef POLYMESH_HH
#define POLYMESH_HH

#pragma once

#include <vector>
#include <string>
#include <unordered_map>

#include "../Math/vectorMath.hh"
#include "../DataBase/field.hh"

namespace Mesh {
    template <int dim>
    class PolyMesh {
    private:
        Field<Lin::Vector<dim>> positions;
        std::vector<std::vector<size_t>> cells;

        std::unordered_map<size_t, std::vector<size_t>> neighbors;
        std::vector<size_t> boundaryNodes;

        std::unordered_map<size_t, std::vector<size_t>> nodeToCellMap;
        std::vector<std::vector<size_t>> cellAdjacency;

    public:
        using Vector = Lin::Vector<dim>;
        using VectorField = Field<Vector>;
        using ScalarField = Field<double>;

        PolyMesh();
        virtual ~PolyMesh() = default;

        void addNode(const Vector& position);
        void setNodePosition(size_t nodeId, const Vector& position);
        size_t addCell(const std::vector<size_t>& nodeIndices);

        const std::vector<Vector>& getNodes() const;
        const std::vector<std::vector<size_t>>& getConnectivityMap() const;

        void computeNeighbors();
        void identifyBoundaryNodes();

        std::vector<size_t> getNeighbors(size_t nodeId) const;
        const std::vector<size_t>& getBoundaryNodes() const;

        void computeConnectivityMap();
        const std::unordered_map<size_t, std::vector<size_t>>& getNodeToCellMap() const;

        void computeCellAdjacency();
        const std::vector<std::vector<size_t>>& getCellAdjacency() const;

        void writeVTK(const std::string& filepath) const;

        // Per-cell scalar measure (area in 2D, volume in 3D) used for the VTK
        // "area" cell-data field. FEMesh answers via Element::computeArea;
        // ALEMesh answers via its own cached polygon area/volume.
        virtual double cellMeasure(size_t cellIndex) const = 0;

        PolyMesh(const PolyMesh& other) = default;
    };
}

#include "polymesh.cc"

#endif // POLYMESH_HH
