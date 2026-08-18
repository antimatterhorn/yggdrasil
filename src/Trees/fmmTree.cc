// Copyright (C) 2026  Cody Raskin

#ifndef FMMTREE_CC
#define FMMTREE_CC

#include <algorithm>
#include <cmath>
#include "fmmTree.hh"

template <int dim>
int FMMTreeNode<dim>::childIndexFor(const Vector& p) const {
    int idx = 0;
    for (int d = 0; d < dim; ++d)
        if (p[d] >= boxCenter[d]) idx |= (1 << d);
    return idx;
}

template <int dim>
typename FMMTreeNode<dim>::Vector FMMTreeNode<dim>::childCenter(int childIdx) const {
    Vector c = boxCenter;
    double q = halfSize * 0.5;
    for (int d = 0; d < dim; ++d)
        c[d] += (childIdx & (1 << d)) ? q : -q;
    return c;
}

template <int dim>
bool FMMTreeNode<dim>::adjacentTo(const FMMTreeNode<dim>& other) const {
    double eps = 1e-9 * std::max(halfSize, other.halfSize);
    for (int d = 0; d < dim; ++d)
        if (std::abs(boxCenter[d] - other.boxCenter[d]) > halfSize + other.halfSize + eps)
            return false;
    return true;
}

template <int dim>
FMMTree<dim>::FMMTree(const Field<Vector>* pos, const Field<double>* mass,
                       int maxSourcesPerLeaf, int maxDepth)
    : positions(pos), masses(mass), maxSourcesPerLeaf(maxSourcesPerLeaf), maxDepth(maxDepth) {}

template <int dim>
void FMMTree<dim>::build() {
    unsigned int n = positions->size();
    if (n == 0) { root.reset(); return; }

    Vector lo = positions->getValue(0);
    Vector hi = lo;
    for (unsigned int i = 1; i < n; ++i) {
        Vector p = positions->getValue(i);
        for (int d = 0; d < dim; ++d) {
            lo[d] = std::min(lo[d], p[d]);
            hi[d] = std::max(hi[d], p[d]);
        }
    }

    Vector center = (lo + hi) * 0.5;
    double halfSize = 1e-12;
    for (int d = 0; d < dim; ++d)
        halfSize = std::max(halfSize, (hi[d] - lo[d]) * 0.5);
    halfSize *= 1.001;

    root = std::make_unique<Node>(center, halfSize, 0);
    root->bodyIndices.resize(n);
    for (unsigned int i = 0; i < n; ++i) root->bodyIndices[i] = i;

    split(root.get());
    subdivideAndClassifySelf(root.get());

    leafOf.assign(n, nullptr);
    assignLeafOf(root.get());
}

template <int dim>
void FMMTree<dim>::assignLeafOf(Node* node) {
    if (node->leaf) {
        for (int idx : node->bodyIndices) leafOf[idx] = node;
        return;
    }
    for (auto& c : node->children)
        if (c) assignLeafOf(c.get());
}

template <int dim>
void FMMTree<dim>::split(Node* node) {
    if ((int)node->bodyIndices.size() <= maxSourcesPerLeaf || node->level >= maxDepth) {
        node->leaf = true;
        return;
    }

    node->leaf = false;
    std::array<std::vector<int>, Node::NumChildren> buckets;
    for (int idx : node->bodyIndices)
        buckets[node->childIndexFor(positions->getValue(idx))].push_back(idx);
    node->bodyIndices.clear();

    for (int c = 0; c < Node::NumChildren; ++c) {
        if (buckets[c].empty()) continue;
        auto child = std::make_unique<Node>(node->childCenter(c), node->halfSize * 0.5, node->level + 1);
        child->parent = node;
        child->bodyIndices = std::move(buckets[c]);
        split(child.get());
        node->children[c] = std::move(child);
    }
}

template <int dim>
void FMMTree<dim>::subdivideAndClassifySelf(Node* node) {
    if (node->leaf) return;

    for (int i = 0; i < Node::NumChildren; ++i) {
        if (!node->children[i]) continue;
        for (int j = i + 1; j < Node::NumChildren; ++j) {
            if (!node->children[j]) continue;
            classifyPair(node->children[i].get(), node->children[j].get());
        }
    }
    for (int i = 0; i < Node::NumChildren; ++i)
        if (node->children[i]) subdivideAndClassifySelf(node->children[i].get());
}

template <int dim>
void FMMTree<dim>::classifyPair(Node* A, Node* B) {
    if (A->adjacentTo(*B)) {
        if (A->leaf && B->leaf) {
            A->uList.push_back(B);
            B->uList.push_back(A);
        } else if (A->leaf) {
            for (int j = 0; j < Node::NumChildren; ++j)
                if (B->children[j]) classifyPair(A, B->children[j].get());
        } else if (B->leaf) {
            for (int i = 0; i < Node::NumChildren; ++i)
                if (A->children[i]) classifyPair(A->children[i].get(), B);
        } else {
            for (int i = 0; i < Node::NumChildren; ++i) {
                if (!A->children[i]) continue;
                for (int j = 0; j < Node::NumChildren; ++j)
                    if (B->children[j]) classifyPair(A->children[i].get(), B->children[j].get());
            }
        }
    } else if (A->level == B->level) {
        A->vList.push_back(B);
        B->vList.push_back(A);
    } else if (A->level < B->level) {
        A->wList.push_back(B);
        B->xList.push_back(A);
    } else {
        B->wList.push_back(A);
        A->xList.push_back(B);
    }
}

template <int dim>
void FMMTree<dim>::upwardPass() {
    if (root) upwardRecurse(root.get());
}

template <int dim>
void FMMTree<dim>::upwardRecurse(Node* node) {
    if (node->leaf) {
        std::vector<Vector> pos;
        std::vector<double> mass;
        pos.reserve(node->bodyIndices.size());
        mass.reserve(node->bodyIndices.size());
        for (int idx : node->bodyIndices) {
            pos.push_back(positions->getValue(idx));
            mass.push_back(masses->getValue(idx));
        }
        node->multipole = CartesianMultipole<dim>::fromSources(pos, mass);
        return;
    }

    std::vector<const CartesianMultipole<dim>*> childMultipoles;
    for (auto& c : node->children) {
        if (!c) continue;
        upwardRecurse(c.get());
        childMultipoles.push_back(&c->multipole);
    }
    node->multipole = CartesianMultipole<dim>::combine(childMultipoles);
}

template <int dim>
void FMMTree<dim>::downwardPass(double G) {
    if (!root) return;
    root->local = LocalExpansion<dim>();
    root->local.center = root->boxCenter;
    accumulateOwnLists(root.get(), G);
    downwardRecurse(root.get(), G);
}

template <int dim>
void FMMTree<dim>::accumulateOwnLists(Node* node, double G) {
    for (auto* src : node->vList)
        node->local += LocalExpansion<dim>::fromMultipole(node->boxCenter, src->multipole, G);
    for (auto* src : node->xList)
        node->local += LocalExpansion<dim>::fromMultipole(node->boxCenter, src->multipole, G);
}

template <int dim>
void FMMTree<dim>::downwardRecurse(Node* node, double G) {
    if (node->leaf) return;
    for (auto& child : node->children) {
        if (!child) continue;
        child->local = node->local.shift(child->boxCenter);
        accumulateOwnLists(child.get(), G);
        downwardRecurse(child.get(), G);
    }
}

template <int dim>
typename FMMTree<dim>::Vector FMMTree<dim>::computeForceOn(int index, double G, double eps2) const {
    Node* leaf = leafOf[index];
    Vector p = positions->getValue(index);

    Vector force = leaf->local.evaluateForce(p);

    for (auto* w : leaf->wList)
        force += LocalExpansion<dim>::fromMultipole(p, w->multipole, G).evaluateForce(p);

    for (int j : leaf->bodyIndices) {
        if (j == index) continue;
        Vector rij = positions->getValue(j) - p;
        force += rij.normal() * (G * masses->getValue(j) / (rij.mag2() + eps2));
    }

    for (auto* u : leaf->uList) {
        for (int j : u->bodyIndices) {
            Vector rij = positions->getValue(j) - p;
            force += rij.normal() * (G * masses->getValue(j) / (rij.mag2() + eps2));
        }
    }

    return force;
}

#endif
