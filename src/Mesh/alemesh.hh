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
        std::vector<double> cellVolumes;

        double computeCellArea(size_t cellIndex) const;
        void computeFaceGeometry(Face<dim>& f) const;

    public:
        using Vector = Lin::Vector<dim>;
        using VectorField = Field<Vector>;
        using ScalarField = Field<double>;

        ALEMesh();

        // Builds the face list from the current cell adjacency (interior
        // faces only -- faces on the outer mesh boundary, where only one
        // cell touches an edge, are not yet constructed; see class comment).
        void computeFaces();

        // Recomputes face normals/areas/centroids and cell volumes from the
        // current node positions. Call after moving nodes; does not change
        // face/cell topology, only geometry.
        void updateGeometry();

        const std::vector<Face<dim>>& getFaces() const;
        double cellVolume(size_t cellIndex) const;
        Vector cellCentroid(size_t cellIndex) const;

        virtual double cellMeasure(size_t cellIndex) const override;

        ALEMesh(const ALEMesh& other) = default;
    };
}

#include "alemesh.cc"

#endif // ALEMESH_HH
