#pragma once

#include "../Math/vectorMath.hh"
#include "../DataBase/field.hh"


namespace Mesh { 

template <int dim>
class VoronoiCell {
    public:
        using Vector = Lin::Vector<dim>;
        Vector generator;
        std::vector<Vector> vertices;
        
        VoronoiCell(Vector generator, std::vector<Vector> vertices)
            : generator(generator),
              vertices(std::move(vertices))
        {}

        Vector getGenerator() const { return generator; }
        const std::vector<Vector>& getVertices() const { return vertices; }
};
  

template <int dim>
class VoronoiMesh {
private:
    using Vector = Lin::Vector<dim>;

    Field<Vector> generators;
    std::vector<VoronoiCell<dim>> cells;



    // Core helpers:
    std::vector<std::pair<Vector, Vector>> getBisectingHyperplanes(int i);
    std::vector<Vector> intersectHalfSpaces(const std::vector<std::pair<Vector, Vector>>& planes);

public:

    VoronoiMesh() = default;

    explicit VoronoiMesh(const Field<Vector>& seedPoints);

    void generateMesh();

    const std::vector<VoronoiCell<dim>>& getCells() const { return cells; }


};
}