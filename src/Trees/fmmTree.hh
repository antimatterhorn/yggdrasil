// Copyright (C) 2026  Cody Raskin

#pragma once

#include <array>
#include <vector>
#include <memory>
#include "../Math/vectorMath.hh"
#include "../DataBase/field.hh"
#include "cartesianMultipole.hh"

// Adaptive box tree with FMM interaction lists (U/V/W/X) and multipole/local expansion evaluation.
template <int dim>
class FMMTreeNode {
public:
    using Vector = Lin::Vector<dim>;
    static constexpr int NumChildren = 1 << dim;

    Vector boxCenter;
    double halfSize;
    int level;
    bool leaf = true;
    FMMTreeNode* parent = nullptr;
    std::vector<int> bodyIndices;
    std::array<std::unique_ptr<FMMTreeNode>, NumChildren> children;

    std::vector<FMMTreeNode*> uList;
    std::vector<FMMTreeNode*> vList;
    std::vector<FMMTreeNode*> wList;
    std::vector<FMMTreeNode*> xList;

    CartesianMultipole<dim> multipole;
    LocalExpansion<dim> local;

    FMMTreeNode(const Vector& boxCenter, double halfSize, int level)
        : boxCenter(boxCenter), halfSize(halfSize), level(level) {}

    int childIndexFor(const Vector& p) const;
    Vector childCenter(int childIdx) const;
    bool adjacentTo(const FMMTreeNode<dim>& other) const;
};

template <int dim>
class FMMTree {
public:
    using Vector = Lin::Vector<dim>;
    using Node = FMMTreeNode<dim>;

    const Field<Vector>* positions;
    const Field<double>* masses;
    int maxSourcesPerLeaf;
    int maxDepth;
    std::unique_ptr<Node> root;
    std::vector<Node*> leafOf;

    FMMTree(const Field<Vector>* pos, const Field<double>* mass,
            int maxSourcesPerLeaf = 16, int maxDepth = 40);

    void build();
    void upwardPass();
    void downwardPass(double G);
    Vector computeForceOn(int index, double G, double eps2) const;

private:
    void split(Node* node);
    void subdivideAndClassifySelf(Node* node);
    void classifyPair(Node* A, Node* B);
    void assignLeafOf(Node* node);
    void upwardRecurse(Node* node);
    void accumulateOwnLists(Node* node, double G);
    void downwardRecurse(Node* node, double G);
};

#include "fmmTree.cc"
