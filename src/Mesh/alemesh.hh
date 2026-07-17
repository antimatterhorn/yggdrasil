// Copyright (C) 2026  Cody Raskin

#ifndef ALEMESH_HH
#define ALEMESH_HH

#pragma once

#include <vector>

#include "../Math/vectorMath.hh"
#include "../DataBase/field.hh"
#include "polymesh.hh"
#include "face.hh"

namespace Mesh {
    template <int dim>
    class ALEMesh : public PolyMesh<dim> {
    private:
        std::vector<Face<dim>> faces;
        std::vector<std::vector<size_t>> facesPerCell;
        std::vector<double> cellVolumes;

        double computeCellArea(size_t cellIndex) const;
        void computeFaceGeometry(Face<dim>& f) const;

    public:
        using Vector = Lin::Vector<dim>;
        using VectorField = Field<Vector>;
        using ScalarField = Field<double>;

        ALEMesh();

        // Builds the face list from the current cell adjacency: one interior
        // Face per pair of adjacent cells, plus one boundary Face (rightCell
        // == Face<dim>::boundaryCell, normal pointing outward) for every cell
        // edge not shared with another cell.
        void computeFaces();

        // Recomputes face normals/areas/centroids and cell volumes from the
        // current node positions. Call after moving nodes; does not change
        // face/cell topology, only geometry.
        void updateGeometry();

        const std::vector<Face<dim>>& getFaces() const;
        // Indices into getFaces() of every face (interior or boundary)
        // touching the given cell. Populated by computeFaces(); lets a hydro
        // package parallelize its derivative loop over cells (each visiting
        // only its own faces) without a scatter/race into shared accumulators.
        const std::vector<size_t>& getFacesForCell(size_t cellIndex) const;
        double cellVolume(size_t cellIndex) const;
        Vector cellCentroid(size_t cellIndex) const;

        virtual double cellMeasure(size_t cellIndex) const override;

        ALEMesh(const ALEMesh& other) = default;
    };
}

#include "alemesh.cc"

#endif // ALEMESH_HH
