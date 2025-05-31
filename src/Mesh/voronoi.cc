#include "voronoi.hh"

namespace Mesh {  

    
template <int dim>
std::vector<std::pair<Lin::Vector<dim>, Lin::Vector<dim>>>
VoronoiMesh<dim>::getBisectingHyperplanes(int i) {
    const auto& pi = generators[i];
    std::vector<std::pair<Lin::Vector<dim>, Lin::Vector<dim>>> planes;

    for (size_t j = 0; j < generators.size(); ++j) {
        if (j == i) continue;
        const auto& pj = generators[j];
        auto midpoint = (pi + pj) * 0.5;
        auto normal = (pj - pi).normal();  // points from pi to pj
        planes.emplace_back(midpoint, -normal);  // inward-facing
    }

    return planes;
}

template <>
std::vector<Lin::Vector<2>>
VoronoiMesh<2>::intersectHalfSpaces(const std::vector<std::pair<Lin::Vector<2>, Lin::Vector<2>>>& planes) {
    using Vector = Lin::Vector<2>;

    // Start with a large square polygon
    std::vector<Vector> polygon = {
        Vector{-1e6, -1e6},
        Vector{+1e6, -1e6},
        Vector{+1e6, +1e6},
        Vector{-1e6, +1e6}
    };

    for (const auto& [p, n] : planes) {
        std::vector<Vector> newPoly;
        for (size_t i = 0; i < polygon.size(); ++i) {
            Vector a = polygon[i];
            Vector b = polygon[(i + 1) % polygon.size()];

            double da = (a - p).dot(n);
            double db = (b - p).dot(n);

            if (da <= 0) newPoly.push_back(a);  // Inside

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
        VoronoiCell cell;
        cell.generator = generators[i];

        auto planes = getBisectingHyperplanes(i);
        cell.vertices = intersectHalfSpaces(planes);

        cells.push_back(std::move(cell));
    }
}

}