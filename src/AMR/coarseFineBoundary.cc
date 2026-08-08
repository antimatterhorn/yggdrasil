// Copyright (C) 2026  Cody Raskin

#pragma once
#include <vector>
#include "../Boundaries/boundary.hh"
#include "../DataBase/nodeList.hh"

namespace AMR {

// Fills a fine patch's ghost rim from its coarse parent by piecewise-constant injection (pairing computed in AMRController.py).
template <int dim>
class CoarseFineBoundary : public Boundary<dim> {
public:
    using Vector       = Lin::Vector<dim>;
    using VectorField  = Field<Vector>;
    using ScalarField  = Field<double>;
    using Complex      = std::complex<double>;
    using ComplexField = Field<Complex>;

    CoarseFineBoundary(NodeList* coarseNodeList,
                        std::vector<int> ghostIds,
                        std::vector<int> coarseIds)
        : coarseNodeList(coarseNodeList),
          ghostIds(std::move(ghostIds)),
          coarseIds(std::move(coarseIds)) {}

    virtual ~CoarseFineBoundary() = default;

    void ApplyScalar(ScalarField* field) override { ApplyThis<double>(field); }
    void ApplyVector(VectorField* field) override { ApplyThis<Vector>(field); }
    void ApplyComplex(ComplexField* field) override { ApplyThis<Complex>(field); }

private:
    NodeList* coarseNodeList;
    std::vector<int> ghostIds;
    std::vector<int> coarseIds;

    template <typename T>
    void ApplyThis(Field<T>* field) {
        Field<T>* coarseField = coarseNodeList->template getFieldOrThrow<T>(field->getName().name());
        #pragma omp parallel for
        for (int k = 0; k < (int)ghostIds.size(); ++k)
            field->setValue(ghostIds[k], coarseField->getValue(coarseIds[k]));
    }
};

}
