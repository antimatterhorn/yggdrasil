// Copyright (C) 2026  Cody Raskin

#pragma once
#include <vector>
#include "../Boundaries/boundary.hh"
#include "../DataBase/nodeList.hh"

namespace AMR {

// Fills a patch's ghost rim from a same-resolution sibling's interior cells (direct copy; pairing computed in AMRController.py).
template <int dim>
class PatchNeighborBoundary : public Boundary<dim> {
public:
    using Vector       = Lin::Vector<dim>;
    using VectorField  = Field<Vector>;
    using ScalarField  = Field<double>;
    using Complex      = std::complex<double>;
    using ComplexField = Field<Complex>;

    PatchNeighborBoundary(NodeList* neighborNodeList,
                           std::vector<int> myGhostIds,
                           std::vector<int> neighborInteriorIds)
        : neighborNodeList(neighborNodeList),
          myGhostIds(std::move(myGhostIds)),
          neighborInteriorIds(std::move(neighborInteriorIds)) {}

    virtual ~PatchNeighborBoundary() = default;

    void ApplyScalar(ScalarField* field) override { ApplyThis<double>(field); }
    void ApplyVector(VectorField* field) override { ApplyThis<Vector>(field); }
    void ApplyComplex(ComplexField* field) override { ApplyThis<Complex>(field); }

private:
    NodeList* neighborNodeList;
    std::vector<int> myGhostIds;
    std::vector<int> neighborInteriorIds;

    template <typename T>
    void ApplyThis(Field<T>* field) {
        Field<T>* neighborField = neighborNodeList->template getFieldOrThrow<T>(field->getName().name());
        #pragma omp parallel for
        for (int k = 0; k < (int)myGhostIds.size(); ++k)
            field->setValue(myGhostIds[k], neighborField->getValue(neighborInteriorIds[k]));
    }
};

}
