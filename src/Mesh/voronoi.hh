// Copyright (C) 2025  Cody Raskin

#pragma once

#include "../Math/vectorMath.hh"
#include "../DataBase/field.hh"
#include "femesh.hh"

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

    double 
    area() const {
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

    Vector centroid() const {
        if constexpr (dim == 2) {
            double cx = 0.0;
            double cy = 0.0;
            double A = 0.0;

            const size_t N = vertices.size();
            if (N < 3) return generator;

            for (size_t i = 0; i < N; ++i) {
                const Vector& p0 = vertices[i];
                const Vector& p1 = vertices[(i + 1) % N];

                double cross = p0.x() * p1.y() - p1.x() * p0.y();
                cx += (p0.x() + p1.x()) * cross;
                cy += (p0.y() + p1.y()) * cross;
                A += cross;
            }

            A *= 0.5;
            if (std::abs(A) < 1e-12) return generator;

            cx /= (6.0 * A);
            cy /= (6.0 * A);
            Vector result{cx, cy};
            return result;
        } else {
            // fallback for unsupported dimensions
            return Vector();  // default-constructed zero vector
        }
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
    FEMesh<dim> generateDualFEMesh() const;

};

}