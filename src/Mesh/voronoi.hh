#pragma once

#include "../Math/vectorMath.hh"

namespace Mesh {  
template <int dim>
class VoronoiMesh {
private:
    using Vector = Lin::Vector<dim>;
    struct VoronoiCell {
        Vector generator;
        std::vector<Vector> vertices; // Convex polytope vertices
        // Optionally, store faces for dim > 2
    };

    std::vector<Vector> generators;
    std::vector<VoronoiCell> cells;



    // Core helpers:
    std::vector<std::pair<Vector, Vector>> getBisectingHyperplanes(int i);
    std::vector<Vector> intersectHalfSpaces(const std::vector<std::pair<Vector, Vector>>& planes);

public:

    VoronoiMesh() = default;

    explicit VoronoiMesh(const std::vector<Vector>& seedPoints)
        : generators(seedPoints) {}

    void generateMesh();

    const std::vector<VoronoiCell>& getCells() const { return cells; }


};
}