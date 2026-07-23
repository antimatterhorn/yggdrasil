// Copyright (C) 2026  Cody Raskin

#ifndef FEMESH_CC
#define FEMESH_CC

#include "femesh.hh"
#include <stdexcept>
#include <set>
#include "../IO/importObj.cc"
#include <iostream>
#include <fstream>
#include <sstream>

namespace Mesh {

template <int dim>
FEMesh<dim>::FEMesh() : PolyMesh<dim>() {}

template <int dim>
void
FEMesh<dim>::addElement(ElementType type, const std::vector<size_t>& nodeIndices) {
    // PolyMesh::addCell validates the node indices and registers the generic
    // connectivity; the Element we construct here mirrors it 1:1 at the same
    // index, so cells[i] (base) and elements[i] (here) always agree.
    this->addCell(nodeIndices);

    auto element = createElement<dim>(type, nodeIndices);
    elements.push_back(element);
    elementTypes.push_back(type);
}

template <int dim>
const std::vector<std::shared_ptr<Element<dim>>>& FEMesh<dim>::getElements() const {
    return elements;
}

template <int dim>
const std::vector<ElementType>& FEMesh<dim>::getElementTypes() const {
    return elementTypes;
}

template <int dim>
std::vector<std::vector<size_t>> FEMesh<dim>::getElementConnectivity() const {
    std::vector<std::vector<size_t>> connectivity;
    for (const auto& elem : elements) {
        connectivity.push_back(elem->nodeIndices());
    }
    return connectivity;
}

template <int dim>
std::vector<std::pair<ElementType, std::vector<size_t>>>
FEMesh<dim>::getElementInfo() const {
    std::vector<std::pair<ElementType, std::vector<size_t>>> result;
    for (size_t i = 0; i < elements.size(); ++i) {
        result.emplace_back(elementTypes[i], elements[i]->nodeIndices());
    }
    return result;
}


void
buildFromObj2dHelper(FEMesh<2>& mesh, const std::string& filepath, const std::string& axes) {
    auto [vertices, faces] = importObj2d(filepath, axes);

    for (const auto& v : vertices.getValues()) {
        mesh.addNode(v);
    }
    int tris=0,quads=0;
    for (const auto& face : faces.getValues()) {
        if (face.size() == 3) {
            tris++;
            mesh.addElement(ElementType::Triangle, {
                static_cast<size_t>(face[0] - 1),
                static_cast<size_t>(face[1] - 1),
                static_cast<size_t>(face[2] - 1)
            });
        } else if (face.size() == 4) {
            quads++;
            mesh.addElement(ElementType::Quad, {
                static_cast<size_t>(face[0] - 1),
                static_cast<size_t>(face[1] - 1),
                static_cast<size_t>(face[2] - 1),
                static_cast<size_t>(face[3] - 1)
            });
        } else {
            throw std::runtime_error("Only triangle and quad elements are supported in buildFromObj.");
        }
    }

    std::cout << "added " << vertices.size() << " verts and " << faces.size() << " faces." << std::endl;
    std::cout << tris << " triangles." << std::endl;
    std::cout << quads << " quads." << std::endl;
}

template <int dim>
void
FEMesh<dim>::buildFromObj(const std::string& filepath, const std::string& axes) {
    std::cout << "loading " << filepath << "..." << std::endl;

    if constexpr (dim == 2) {
        buildFromObj2dHelper(*this, filepath, axes);
    } else {
        throw std::runtime_error("buildFromObj is only implemented for 2D FEMesh");
    }

    this->computeConnectivityMap();
}

template <int dim>
double
FEMesh<dim>::cellMeasure(size_t cellIndex) const {
    return elements[cellIndex]->computeArea(this->getNodes());
}

// Explicit instantiations
template class FEMesh<2>;
template class FEMesh<3>;

} // namespace Mesh

#endif // FEMESH_CC
