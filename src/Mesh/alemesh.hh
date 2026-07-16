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
    // Finite-volume mesh substrate for ALE hydro: PolyMesh's generic node/cell
    // connectivity plus a Face list (normal, area, centroid, adjacent cells) and
    // per-cell volumes, both recomputed from current node positions on demand
    // via updateGeometry(). Cells are added as plain node-index polygons via
    // the inherited PolyMesh::addCell -- no typed Element/basis-function
    // machinery, since finite-volume flux divergence has no use for it.
    //
    // 2D only for now: face construction assumes a 2-node (edge) shared
    // boundary between adjacent cells, matching Element's current 2D-only
    // scope elsewhere in this module. A hydro solver that actually walks these
    // faces, plus mesh-velocity/remap support, is separate follow-up work.
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
