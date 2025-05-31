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
        if (j == i) continue;
        const auto& pj = generators.getValue(j);
        auto midpoint = (pi + pj) * 0.5;
        auto normal = (pj - pi).normal();  // points from pi to pj
        planes.emplace_back(midpoint, -normal);  // inward-facing
    }

    return planes;
}

template <>
std::vector<Lin::Vector<3>>
VoronoiMesh<3>::intersectHalfSpaces(const std::vector<std::pair<Lin::Vector<3>, Lin::Vector<3>>>& planes) {
    throw std::runtime_error("3D Voronoi not implemented yet");
}

template <>
std::vector<Lin::Vector<2>>
VoronoiMesh<2>::intersectHalfSpaces(const std::vector<std::pair<Lin::Vector<2>, Lin::Vector<2>>>& planes) {
    using Vector = Lin::Vector<2>;

    // Compute bounding box of generator points
    Vector bmin = generators.getValue(0);
    Vector bmax = generators.getValue(0);
    for (size_t i = 1; i < generators.size(); ++i) {
        const Vector& p = generators.getValue(i);
        for (int d = 0; d < 2; ++d) {
            bmin[d] = std::min(bmin[d], p[d]);
            bmax[d] = std::max(bmax[d], p[d]);
        }
    }

    // Pad by 50% margin
    Vector size = bmax - bmin;
    Vector pad = size * 0.5;
    bmin -= pad;
    bmax += pad;

    std::vector<Vector> polygon = {
        Vector{bmin[0], bmin[1]},
        Vector{bmax[0], bmin[1]},
        Vector{bmax[0], bmax[1]},
        Vector{bmin[0], bmax[1]}
    };

    for (const auto& [p, n] : planes) {
        std::vector<Vector> newPoly;
        for (size_t i = 0; i < polygon.size(); ++i) {
            Vector a = polygon[i];
            Vector b = polygon[(i + 1) % polygon.size()];

            double da = (a - p).dot(n);
            double db = (b - p).dot(n);

            if (da <= 0) newPoly.push_back(a);

            if ((da <= 0 && db > 0) || (da > 0 && db <= 0)) {
                double t = da / (da - db);
                Vector intersection = a + (b - a) * t;
                newPoly.push_back(intersection);
            }
        }
        polygon = std::move(newPoly);
        if (polygon.empty()) break;
    }

    return polygon;
}


template <int dim>
void VoronoiMesh<dim>::generateMesh() {
    cells.clear();
    for (size_t i = 0; i < generators.size(); ++i) {
        auto planes = getBisectingHyperplanes(i);
        VoronoiCell<dim> cell(generators.getValue(i),intersectHalfSpaces(planes));

        cells.push_back(std::move(cell));
    }
}

}

template class Mesh::VoronoiCell<2>;
template class Mesh::VoronoiCell<3>;