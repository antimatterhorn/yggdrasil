// Copyright (C) 2026  Cody Raskin

#ifndef ALEMESH_CC
#define ALEMESH_CC

#include "alemesh.hh"
#include <set>
#include <utility>
#include <algorithm>

namespace Mesh {

template <int dim>
ALEMesh<dim>::ALEMesh() : PolyMesh<dim>() {}

template <int dim>
typename ALEMesh<dim>::Vector
ALEMesh<dim>::cellCentroid(size_t cellIndex) const {
    const auto& nodeIds = this->getConnectivityMap()[cellIndex];
    const auto& nodes = this->getNodes();
    Vector c = Lin::Vector<dim>::zero();
    for (size_t id : nodeIds) c = c + nodes[id];
    return c * (1.0 / static_cast<double>(nodeIds.size()));
}

template <int dim>
double
ALEMesh<dim>::computeCellArea(size_t cellIndex) const {
    // Shoelace formula over the cell's current node positions.
    if constexpr (dim == 2) {
        const auto& nodeIds = this->getConnectivityMap()[cellIndex];
        const auto& nodes = this->getNodes();
        double A = 0.0;
        size_t N = nodeIds.size();
        for (size_t i = 0; i < N; ++i) {
            const Vector& p1 = nodes[nodeIds[i]];
            const Vector& p2 = nodes[nodeIds[(i + 1) % N]];
            A += p1.x() * p2.y() - p2.x() * p1.y();
        }
        return 0.5 * std::abs(A);
    } else {
        return 0.0;
    }
}

template <int dim>
void
ALEMesh<dim>::computeFaceGeometry(Face<dim>& f) const {
    if constexpr (dim == 2) {
        if (f.nodeIndices.size() != 2) return;
        const auto& nodes = this->getNodes();
        Vector p0 = nodes[f.nodeIndices[0]];
        Vector p1 = nodes[f.nodeIndices[1]];
        Vector edge = p1 - p0;

        f.area = edge.magnitude();
        f.centroid = (p0 + p1) * 0.5;

        Vector n = Lin::Vector<dim>::zero();
        n.setX(-edge.y());
        n.setY(edge.x());
        double nm = n.magnitude();
        if (nm > 1e-14) n = n * (1.0 / nm);

        // Orient the normal to point from leftCell toward rightCell: if it
        // currently points back toward leftCell's own centroid instead of
        // away from it, flip it. A consistent orientation is required for
        // any future flux divergence built on top of this.
        if (f.leftCell != Face<dim>::boundaryCell) {
            Vector lc = cellCentroid(f.leftCell);
            if (n.dot(f.centroid - lc) < 0.0) n = n * -1.0;
        }
        f.normal = n;
    }
}

template <int dim>
void
ALEMesh<dim>::computeFaces() {
    this->computeConnectivityMap();
    this->computeCellAdjacency();

    const auto& adjacency = this->getCellAdjacency();
    const auto& cells = this->getConnectivityMap();

    faces.clear();
    std::set<std::pair<size_t, size_t>> seen;

    for (size_t ci = 0; ci < adjacency.size(); ++ci) {
        for (size_t cj : adjacency[ci]) {
            auto key = std::make_pair(std::min(ci, cj), std::max(ci, cj));
            if (!seen.insert(key).second) continue;

            std::vector<size_t> shared;
            for (size_t a : cells[ci])
                for (size_t b : cells[cj])
                    if (a == b) shared.push_back(a);

            Face<dim> f;
            f.nodeIndices = shared;
            f.leftCell = ci;
            f.rightCell = cj;
            computeFaceGeometry(f);
            faces.push_back(f);
        }
    }

    cellVolumes.resize(cells.size());
    for (size_t ci = 0; ci < cells.size(); ++ci)
        cellVolumes[ci] = computeCellArea(ci);
}

template <int dim>
void
ALEMesh<dim>::updateGeometry() {
    const auto& cells = this->getConnectivityMap();
    cellVolumes.resize(cells.size());
    for (size_t ci = 0; ci < cells.size(); ++ci)
        cellVolumes[ci] = computeCellArea(ci);

    for (auto& f : faces)
        computeFaceGeometry(f);
}

template <int dim>
const std::vector<Face<dim>>& ALEMesh<dim>::getFaces() const {
    return faces;
}

template <int dim>
double
ALEMesh<dim>::cellVolume(size_t cellIndex) const {
    return cellMeasure(cellIndex);
}

template <int dim>
double
ALEMesh<dim>::cellMeasure(size_t cellIndex) const {
    if (cellIndex < cellVolumes.size()) return cellVolumes[cellIndex];
    return computeCellArea(cellIndex);
}

// Explicit instantiation: 2D only for now (see class comment in alemesh.hh).
template class ALEMesh<2>;

} // namespace Mesh

#endif // ALEMESH_CC
