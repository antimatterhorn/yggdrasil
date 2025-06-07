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

        double area() const {
        double A = 0.0;
        size_t N = vertices.size();
        if (dim == 2) {
            for (size_t i = 0; i < N; ++i) {
                const Vector& p1 = vertices[i];
                const Vector& p2 = vertices[(i + 1) % N];
                A += p1.x() * p2.y() - p2.x() * p1.y();
            }
        }
        return 0.5 * std::abs(A);
    }
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