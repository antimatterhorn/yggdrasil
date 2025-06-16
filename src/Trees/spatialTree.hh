#pragma once

#include <array>
#include <vector>
#include <memory>
#include "../Math/vectorMath.hh" 
#include "../DataBase/field.hh" 

template<int dim>
class SpatialTreeNode {
public:
    using Vector = Lin::Vector<dim>;
    static constexpr int NumChildren = 1 << dim;

    Vector center;      // Center of this node
    double halfSize;    // Half the size of this node's region
    double mass = 0.0;
    Vector com = Vector::zero();  // Center of mass
    int bodyIndex = -1;           // Index in global array; -1 if internal
    std::array<std::unique_ptr<SpatialTreeNode>, NumChildren> children;

    SpatialTreeNode(Vector center, double halfSize)
        : center(center), halfSize(halfSize) {}

    bool contains(const Vector& p) const {
        for (int d = 0; d < dim; ++d)
            if (std::abs(p[d] - center[d]) > halfSize) return false;
        return true;
    }

    int getChildIndex(const Vector& p) const {
        int idx = 0;
        for (int d = 0; d < dim; ++d)
            if (p[d] >= center[d]) idx |= (1 << d);
        return idx;
    }

    void insert(int index, const Field<Vector>& positions, const Field<double>& masses);

    void computeForce(int index,
                    const Vector& position,
                    double theta,
                    double G,
                    double eps2,
                    Vector& acc,
                    const Field<double>& masses) const;
};

template<int dim>
class SpatialTree {
public:
    using Vector = Lin::Vector<dim>;

    const Field<Vector>* positions = nullptr;
    const Field<double>* masses = nullptr;
    std::unique_ptr<SpatialTreeNode<dim>> root;

    SpatialTree(const Field<Vector>* pos,
                const Field<double>* mass);

    void build();

    Vector computeForceOn(int index, double theta, double G, double eps2) const;
};

#include "spatialTree.cc"

