// Copyright (C) 2025  Cody Raskin

#pragma once

#include "constraint.hh"

// Base class for Grid Boundaries
template <int dim>
class MotionConstraint : public Constraint<dim> {
protected:
    Lin::Vector<dim> direction;

public:
    using Vector=Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;
    
    MotionConstraint(NodeList* nodeList, std::vector<int> nodeIndices, Vector direction) : 
        Constraint<dim>(nodeList, nodeIndices), direction(direction) { }

    virtual ~MotionConstraint() {}

    virtual void
    ApplyBoundaries(State<dim>* state, NodeList* nodeList) override {
        State<dim> copyState = this->copyState;
        std::vector<int> nodeIndices = this->nodeIndices;
        VectorField* position = state->template getField<Vector>("position");
        VectorField* copyPos  = copyState.template getField<Vector>("position");
        VectorField* velocity = state->template getField<Vector>("velocity");
        VectorField* copyVel  = copyState.template getField<Vector>("velocity");
        
        for (int i = 0; i < nodeIndices.size(); ++i) {
            int index = nodeIndices[i];
            Vector newPos = position->getValue(index);
            Vector newVel = velocity->getValue(index);
            for (int d=0; d<dim; d++) {
                newPos[d] = (direction[d] ? copyPos->getValue(i)[d] : position->getValue(index)[d]);
                newVel[d] = (direction[d] ? 0.0 : velocity->getValue(index)[d]);
            }
            position->setValue(index, newPos);
            velocity->setValue(index, newVel);
        }
    }
};