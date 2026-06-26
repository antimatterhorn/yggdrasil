#include "spatialTree.hh"
#include <cmath>
#include <cassert>

//=======================//
// SpatialTreeNode<dim> //
//=======================//

template<int dim>
void SpatialTreeNode<dim>::insert(int index,
                                  const Field<Vector>& positions,
                                  const Field<double>& masses) {
    const Vector& pos = positions.getValue(index);
    const double m = masses.getValue(index);

    // Update center of mass and total mass
    if (mass == 0.0) {
        com = pos;
        mass = m;
    } else {
        com = (com * mass + pos * m) / (mass + m);
        mass += m;
    }

    if (bodyIndex == -1 && !children[0]) {
        // First real insertion
        bodyIndex = index;
        return;
    }

    if (!children[0]) {
        // Subdivide and move existing body
        for (int i = 0; i < NumChildren; ++i) {
            Vector childCenter = center;
            for (int d = 0; d < dim; ++d) {
                childCenter[d] += ((i >> d) & 1 ? 0.5 : -0.5) * halfSize;
            }
            children[i] = std::make_unique<SpatialTreeNode<dim>>(childCenter, halfSize * 0.5);
        }

        if (bodyIndex != -1) {
            int oldChild = getChildIndex(positions.getValue(bodyIndex));
            children[oldChild]->insert(bodyIndex, positions, masses);
            bodyIndex = -1;
        }
    }

    int childIdx = getChildIndex(pos);
    children[childIdx]->insert(index, positions, masses);
}

template<int dim>
void SpatialTreeNode<dim>::computeForce(int index,
                                        const Vector& position,
                                        double theta,
                                        double G,
                                        double eps2,
                                        Vector& acc,
                                        const Field<double>& masses) const {
    if (mass == 0.0 || (bodyIndex == index && !children[0])) return;

    Vector disp = com - position;
    double dist2 = disp.mag2() + eps2;
    double dist = std::sqrt(dist2);

    if (!children[0] || (2.0 * halfSize / dist < theta)) {
        acc += G * mass / (dist2) * (disp / dist);
    } else {
        for (const auto& child : children) {
            if (child) {
                child->computeForce(index, position, theta, G, eps2, acc, masses);
            }
        }
    }
}

//==================//
// SpatialTree<dim> //
//==================//

template<int dim>
SpatialTree<dim>::SpatialTree(const Field<Vector>* pos,
                               const Field<double>* mass)
    : positions(pos), masses(mass) {}

template<int dim>
void SpatialTree<dim>::build() {
    assert(positions && masses);
    assert(positions->size() == masses->size());
    const int N = positions->size();

    // Compute bounding box
    Vector center = Vector::zero();
    for (int i = 0; i < N; ++i) center += positions->getValue(i);
    center = center/N;

    double maxR = 0.0;
    for (int i = 0; i < N; ++i)
        maxR = std::max(maxR, (positions->getValue(i) - center).magnitude());

    root = std::make_unique<SpatialTreeNode<dim>>(center, maxR * 1.1);

    for (int i = 0; i < N; ++i)
        root->insert(i, *positions, *masses);
}

template<int dim>
typename SpatialTree<dim>::Vector
SpatialTree<dim>::computeForceOn(int index,
                                 double theta,
                                 double G,
                                 double eps2) const {
    Vector acc = Vector::zero();
    root->computeForce(index, positions->getValue(index), theta, G, eps2, acc, *masses);
    return acc;
}

//==================//
// Explicit instantiation
//==================//
template class SpatialTreeNode<2>;
template class SpatialTreeNode<3>;
template class SpatialTree<2>;
template class SpatialTree<3>;
