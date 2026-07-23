// Copyright (C) 2026  Cody Raskin

#ifndef FEMESH_HH
#define FEMESH_HH

#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <array>

#include "../Math/vectorMath.hh"
#include "../DataBase/field.hh"
#include "../DataBase/nodeList.hh"
#include "polymesh.hh"
#include "element.hh"

namespace Mesh {
    template <int dim>
    class FEMesh : public PolyMesh<dim> {
    private:
        std::vector<std::shared_ptr<Element<dim>>> elements;
        std::vector<ElementType> elementTypes;

    public:
        using Vector = Lin::Vector<dim>;
        using VectorField = Field<Vector>;
        using ScalarField = Field<double>;

        FEMesh();

        void addElement(ElementType type, const std::vector<size_t>& nodeIndices);

        const std::vector<std::shared_ptr<Element<dim>>>& getElements() const;
        const std::vector<ElementType>& getElementTypes() const;

        std::vector<std::vector<size_t>> getElementConnectivity() const;
        std::vector<std::pair<ElementType, std::vector<size_t>>> getElementInfo() const;

        void buildFromObj(const std::string& filepath, const std::string& axes);

        virtual double cellMeasure(size_t cellIndex) const override;

        FEMesh(const FEMesh& other) = default;
    };
}

#include "femesh.cc"

#endif // FEMESH_HH
