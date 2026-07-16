// Copyright (C) 2026  Cody Raskin

#ifndef POLYMESH_CC
#define POLYMESH_CC

#include "polymesh.hh"
#include <stdexcept>
#include <set>
#include <fstream>
#include <iostream>

namespace Mesh {

template <int dim>
PolyMesh<dim>::PolyMesh() : positions("positions") {}

template <int dim>
void
PolyMesh<dim>::addNode(const Vector& position) {
    positions.addValue(position);
}

template <int dim>
void
PolyMesh<dim>::setNodePosition(size_t nodeId, const Vector& position) {
    positions.setValue(static_cast<unsigned int>(nodeId), position);
}

template <int dim>
size_t
PolyMesh<dim>::addCell(const std::vector<size_t>& nodeIndices) {
    for (size_t idx : nodeIndices) {
        if (idx >= positions.size()) {
            throw std::out_of_range("Node index out of bounds in addCell()");
        }
    }
    cells.push_back(nodeIndices);
    return cells.size() - 1;
}

template <int dim>
const std::vector<typename PolyMesh<dim>::Vector>& PolyMesh<dim>::getNodes() const {
    return positions.getValues();
}

template <int dim>
const std::vector<std::vector<size_t>>& PolyMesh<dim>::getConnectivityMap() const {
    return cells;
}

template <int dim>
void
PolyMesh<dim>::computeNeighbors() {
    neighbors.clear();

    for (const auto& nodeIds : cells) {
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            size_t ni = nodeIds[i];
            for (size_t j = 0; j < nodeIds.size(); ++j) {
                if (i != j) {
                    neighbors[ni].push_back(nodeIds[j]);
                }
            }
        }
    }

    // Remove duplicate neighbors
    for (auto& pair : neighbors) {
        std::set<size_t> unique(pair.second.begin(), pair.second.end());
        pair.second.assign(unique.begin(), unique.end());
    }
}

template <int dim>
std::vector<size_t>
PolyMesh<dim>::getNeighbors(size_t nodeId) const {
    auto it = neighbors.find(nodeId);
    if (it != neighbors.end()) {
        return it->second;
    } else {
        return {};
    }
}

template <int dim>
void
PolyMesh<dim>::identifyBoundaryNodes() {
    // Placeholder: could later be based on cell adjacency, tags, or external criteria
    boundaryNodes.clear();
    for (size_t i = 0; i < positions.size(); ++i) {
        if (neighbors.find(i) == neighbors.end() || neighbors.at(i).size() < 2) {
            boundaryNodes.push_back(i);
        }
    }
}

template <int dim>
const std::vector<size_t>& PolyMesh<dim>::getBoundaryNodes() const {
    return boundaryNodes;
}

template <int dim>
void PolyMesh<dim>::computeConnectivityMap() {
    nodeToCellMap.clear();
    for (size_t ci = 0; ci < cells.size(); ++ci) {
        for (size_t nodeId : cells[ci]) {
            nodeToCellMap[nodeId].push_back(ci);
        }
    }
}

template <int dim>
const std::unordered_map<size_t, std::vector<size_t>>& PolyMesh<dim>::getNodeToCellMap() const {
    return nodeToCellMap;
}

template <int dim>
void PolyMesh<dim>::computeCellAdjacency() {
    // Requires computeConnectivityMap() to have been called first so
    // nodeToCellMap is current.
    cellAdjacency.assign(cells.size(), {});
    for (size_t ci = 0; ci < cells.size(); ++ci) {
        std::unordered_map<size_t, int> shareCount;
        for (size_t nodeId : cells[ci]) {
            auto it = nodeToCellMap.find(nodeId);
            if (it == nodeToCellMap.end()) continue;
            for (size_t cj : it->second) {
                if (cj != ci) shareCount[cj]++;
            }
        }
        // Two cells are adjacent (share a face) if they share at least `dim`
        // node indices -- an edge in 2D, a triangular facet in 3D.
        for (const auto& entry : shareCount) {
            if (entry.second >= dim) {
                cellAdjacency[ci].push_back(entry.first);
            }
        }
    }
}

template <int dim>
const std::vector<std::vector<size_t>>& PolyMesh<dim>::getCellAdjacency() const {
    return cellAdjacency;
}

template <int dim>
void
PolyMesh<dim>::writeVTK(const std::string& filepath) const {
    std::ofstream file(filepath);
    if (!file) throw std::runtime_error("Failed to open file for writing: " + filepath);
    std::cout << "writing to " << filepath << "..." << std::endl;
    file << "# vtk DataFile Version 3.0\n";
    file << "Mesh Export\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    const auto& nodes = positions.getValues();
    file << "POINTS " << nodes.size() << " double\n";
    for (const auto& node : nodes) {
        for (int i = 0; i < dim; ++i) file << node[i] << " ";
        if (dim == 2) file << "0.0";  // VTK expects 3D coords
        file << "\n";
    }

    size_t totalIndices = 0;
    for (const auto& c : cells) totalIndices += c.size();
    file << "CELLS " << cells.size() << " " << (cells.size() + totalIndices) << "\n";
    for (const auto& c : cells) {
        file << c.size();
        for (auto id : c) file << " " << id;
        file << "\n";
    }

    file << "CELL_TYPES " << cells.size() << "\n";
    for (const auto& c : cells) {
        switch (c.size()) {
            case 3: file << "5\n"; break;  // VTK_TRIANGLE
            case 4: file << "9\n"; break;  // VTK_QUAD
            default: file << "0\n"; break; // unknown
        }
    }

    file << "\nPOINT_DATA " << nodes.size() << "\n";
    file << "SCALARS x double 1\nLOOKUP_TABLE default\n";
    for (const auto& node : nodes) file << node[0] << "\n";
    file << "\nSCALARS y double 1\nLOOKUP_TABLE default\n";
    for (const auto& node : nodes) file << ((dim > 1) ? node[1] : 0.0) << "\n";

    file << "\nCELL_DATA " << cells.size() << "\n";
    file << "SCALARS area double 1\nLOOKUP_TABLE default\n";
    for (size_t ci = 0; ci < cells.size(); ++ci) file << cellMeasure(ci) << "\n";
}

// Explicit instantiations
template class PolyMesh<2>;
template class PolyMesh<3>;

} // namespace Mesh

#endif // POLYMESH_CC
