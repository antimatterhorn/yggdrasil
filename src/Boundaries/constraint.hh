// Copyright (C) 2025  Cody Raskin

#pragma once

#include "boundary.hh"
#include "../DataBase/nodeList.hh"
#include "../DataBase/field.hh"
#include <vector>

// Base class for Grid Boundaries
template <int dim>
class Constraint : public Boundary<dim> {
protected:
    NodeList* nodeList;
    std::vector<int> nodeIndices;
    State<dim> copyState;


public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;
    
    Constraint(NodeList* nodeList, std::vector<int> nodeIndices) : 
        nodeList(nodeList), nodeIndices(nodeIndices), copyState(nodeIndices.size()) {
    }

    virtual ~Constraint() {}

    virtual void ZeroTimeInitialize() override {
        std::vector<std::string> names = nodeList->fieldNames();
        
        for (int i = 0; i < names.size(); ++i) {
            const std::string& name = names[i];
            FieldBase* FieldToCopy = nodeList->getFieldByIndex(i);
            
            if (auto* resultDouble = dynamic_cast<Field<double>*>(FieldToCopy)) {
                copyState.template insertField<double>(name);
                ScalarField* copy = copyState.template getField<double>(name);
                for (int j = 0; j < nodeIndices.size(); ++j) {
                    copy->setValue(j, resultDouble->getValue(nodeIndices[j]));
                }
            } 
            else if (auto* resultVector = dynamic_cast<Field<Lin::Vector<dim>>*>(FieldToCopy)) {
                copyState.template insertField<Vector>(name);
                VectorField* copy = copyState.template getField<Vector>(name);
                for (int j = 0; j < nodeIndices.size(); ++j) {
                    copy->setValue(j, resultVector->getValue(nodeIndices[j]));
                }
            }
        }
    }


};