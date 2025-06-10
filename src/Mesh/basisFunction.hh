#pragma once

#include "../Math/vectorMath.hh"
#include <cstddef>
#include <stdexcept>

namespace Mesh {

class BasisFunction {
public:
    using Vector2 = Lin::Vector<2>;

    virtual ~BasisFunction() = default;

    // Evaluate basis function i at reference point (xi, eta)
    virtual double evaluate(size_t i, const Vector2& xi) const = 0;

    // Evaluate gradient in reference coordinates
    virtual Vector2 gradient(size_t i, const Vector2& xi) const = 0;

    virtual size_t numNodes() const = 0;
};

class TriangleBasisFunction : public BasisFunction {
public:
    double evaluate(size_t i, const Vector2& xi) const override {
        double x = xi[0], y = xi[1];
        switch (i) {
            case 0: return 1.0 - x - y;
            case 1: return x;
            case 2: return y;
            default: throw std::out_of_range("Triangle basis index");
        }
    }

    Vector2 gradient(size_t i, const Vector2&) const override {
        switch (i) {
            case 0: return {-1.0, -1.0};
            case 1: return {1.0,  0.0};
            case 2: return {0.0,  1.0};
            default: throw std::out_of_range("Triangle basis grad index");
        }
    }

    size_t numNodes() const override { return 3; }
};

class QuadBasisFunction : public BasisFunction {
public:
    // Node reference coordinates: (-1,-1), (1,-1), (1,1), (-1,1)
    static constexpr double nodeXi[4][2] = {
        {-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}
    };

    double evaluate(size_t i, const Vector2& xi) const override {
        if (i >= 4) throw std::out_of_range("Quad basis index");
        return 0.25 * (1 + xi[0] * nodeXi[i][0]) * (1 + xi[1] * nodeXi[i][1]);
    }

    Vector2 gradient(size_t i, const Vector2& xi) const override {
        if (i >= 4) throw std::out_of_range("Quad basis grad index");

        double dN_dxi = 0.25 * nodeXi[i][0] * (1 + xi[1] * nodeXi[i][1]);
        double dN_deta = 0.25 * (1 + xi[0] * nodeXi[i][0]) * nodeXi[i][1];
        return {dN_dxi, dN_deta};
    }

    size_t numNodes() const override { return 4; }
};

} // namespace Mesh
