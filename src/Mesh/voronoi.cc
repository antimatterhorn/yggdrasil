// Copyright (C) 2025  Cody Raskin

#include "voronoi.hh"
#include <stdexcept>

namespace Mesh {  

template <int dim>
VoronoiMesh<dim>::VoronoiMesh(const Field<Lin::Vector<dim>>& seedPoints)
        : generators(seedPoints) {
            generateMesh();
}

template <int dim>
std::vector<std::pair<Lin::Vector<dim>, Lin::Vector<dim>>>
VoronoiMesh<dim>::getBisectingHyperplanes(int i) {
    const auto& pi = generators.getValue(i);
    std::vector<std::pair<Lin::Vector<dim>, Lin::Vector<dim>>> planes;

    for (size_t j = 0; j < generators.size(); ++j) {
        if (j == (size_t)i) continue;
        const auto& pj = generators.getValue(j);
        auto midpoint = (pi + pj) * 0.5;
        auto normal = (pj - pi).normal();
        planes.emplace_back(midpoint, -normal);  // inward facing
    }

    // Nearest neighbors clip first; distant planes become no-ops early
    std::sort(planes.begin(), planes.end(), [&](const auto& a, const auto& b) {
        return (a.first - pi).mag2() < (b.first - pi).mag2();
    });

    return planes;
}

template <>
std::vector<Lin::Vector<3>>
VoronoiMesh<3>::intersectHalfSpaces(const std::vector<std::pair<Lin::Vector<3>, Lin::Vector<3>>>& planes, const Lin::Vector<3>& bmin, const Lin::Vector<3>& bmax) {
    throw std::runtime_error("3D Voronoi not implemented yet");
}

template <>
std::vector<Lin::Vector<2>>
VoronoiMesh<2>::intersectHalfSpaces(const std::vector<std::pair<Lin::Vector<2>, Lin::Vector<2>>>& planes, const Lin::Vector<2>& bmin, const Lin::Vector<2>& bmax) {
    using Vector = Lin::Vector<2>;

    std::vector<Vector> polygon = {
        Vector{bmin[0], bmin[1]},
        Vector{bmax[0], bmin[1]},
        Vector{bmax[0], bmax[1]},
        Vector{bmin[0], bmax[1]}
    };

    int consecutive_no_clip = 0;
    for (const auto& [p, n] : planes) {
        std::vector<Vector> newPoly;
        bool clipped = false;
        for (size_t i = 0; i < polygon.size(); ++i) {
            Vector a = polygon[i];
            Vector b = polygon[(i + 1) % polygon.size()];

            double da = (a - p).dot(n);
            double db = (b - p).dot(n);

            if (da >= 0) newPoly.push_back(a);
            else clipped = true;

            if ((da >= 0 && db < 0) || (da < 0 && db >= 0)) {
                double t = da / (da - db);
                newPoly.push_back(a + (b - a) * t);
                clipped = true;
            }
        }
        polygon = std::move(newPoly);
        if (polygon.empty()) break;

        // With planes sorted by distance, once several consecutive far planes
        // leave the polygon unchanged the cell is fully defined.
        if (!clipped) {
            if (++consecutive_no_clip >= 5) break;
        } else {
            consecutive_no_clip = 0;
        }
    }

    return polygon;
}


template <int dim>
void VoronoiMesh<dim>::generateMesh() {
    cells.clear();

    // Compute bounding box once; pass it into every intersectHalfSpaces call
    // instead of recomputing it O(N) times inside that function.
    Vector bmin = generators.getValue(0);
    Vector bmax = generators.getValue(0);
    for (size_t i = 1; i < generators.size(); ++i) {
        const Vector& p = generators.getValue(i);
        for (int d = 0; d < dim; ++d) {
            bmin[d] = std::min(bmin[d], p[d]);
            bmax[d] = std::max(bmax[d], p[d]);
        }
    }
    Vector pad = (bmax - bmin) * 0.5;
    bmin -= pad;
    bmax += pad;

    for (size_t i = 0; i < generators.size(); ++i) {
        auto planes = getBisectingHyperplanes(i);
        cells.push_back(VoronoiCell<dim>(generators.getValue(i), intersectHalfSpaces(planes, bmin, bmax)));
    }
}

template <int dim>
FEMesh<dim> VoronoiMesh<dim>::generateDualFEMesh() const {
    if constexpr (dim != 2) {
        throw std::runtime_error("generateDualFEMesh is only implemented for dim = 2");
    } else {
        FEMesh<2> dualMesh;
        const size_t n = generators.size();

        // Add all generator nodes to the dual FEMesh
        for (size_t i = 0; i < n; ++i) {
            dualMesh.addNode(generators[i]);
        }

        // Use quantized (rounded) coordinates to avoid floating-point mismatch
        auto quantize = [](double x) {
            return static_cast<int64_t>(std::round(x * 1e6));  // Adjust precision as needed
        };

        std::map<std::pair<int64_t, int64_t>, std::vector<size_t>> vertexToGenerators;

        // Build a mapping from each Voronoi vertex to all generators (cell indices) that share it
        for (size_t i = 0; i < cells.size(); ++i) {
            for (const auto& v : cells[i].getVertices()) {
                std::pair<int64_t, int64_t> key = {quantize(v.x()), quantize(v.y())};
                vertexToGenerators[key].push_back(i);
            }
        }

        // For each shared Voronoi vertex, build triangles between its associated generators
        for (const auto& [key, genIndices] : vertexToGenerators) {
            if (genIndices.size() >= 3) {
                // Build triangle fan from first generator
                for (size_t i = 1; i + 1 < genIndices.size(); ++i) {
                    std::vector<size_t> tri = {
                        genIndices[0],
                        genIndices[i],
                        genIndices[i + 1]
                    };
                    dualMesh.addElement(ElementType::Triangle, tri);
                }
            }
        }

        std::cout << "Returning FEMesh with " << dualMesh.getNodes().size()
                  << " nodes and " << dualMesh.getElements().size() << " elements\n";

        return dualMesh;
    }
}


}

template class Mesh::VoronoiCell<2>;
template class Mesh::VoronoiCell<3>;

