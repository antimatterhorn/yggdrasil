// Copyright (C) 2026  Cody Raskin

#pragma once

#include "../Math/vectorMath.hh"

namespace Mesh {

    // A face shared between two adjacent ALEMesh cells: the node indices that
    // define it, plus its computed geometry. Normal is oriented to point from
    // leftCell toward rightCell (see ALEMesh::computeFaceGeometry).
    template <int dim>
    struct Face {
        using Vector = Lin::Vector<dim>;

        static constexpr size_t boundaryCell = static_cast<size_t>(-1);

        std::vector<size_t> nodeIndices;
        Vector normal;
        double area = 0.0;
        Vector centroid;

        size_t leftCell = boundaryCell;
        size_t rightCell = boundaryCell;

        bool isBoundary() const { return leftCell == boundaryCell || rightCell == boundaryCell; }

        const std::vector<size_t>& getNodeIndices() const { return nodeIndices; }
        Vector getNormal() const { return normal; }
        double getArea() const { return area; }
        Vector getCentroid() const { return centroid; }
        size_t getLeftCell() const { return leftCell; }
        size_t getRightCell() const { return rightCell; }
    };

}
