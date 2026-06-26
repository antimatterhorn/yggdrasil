// Copyright (C) 2025  Cody Raskin

#ifndef PHYSICS_HH
#define PHYSICS_HH

#include "../Math/vectorMath.hh"
#include "../Type/physicalConstants.hh"
#include "../State/state.hh"
#include "../Boundaries/boundary.hh"

template <int dim>
class Boundary; // forward declaration

template <int dim>
class Physics {
protected:
    NodeList* nodeList;
    PhysicalConstants& constants;
    State<dim> state;
    double lastDt;
    std::vector<Boundary<dim>*> boundaries;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    Physics(NodeList* nodeList, PhysicalConstants& constants) : 
        nodeList(nodeList), 
        constants(constants),
        state(nodeList->size()) {
        VerifyFields(nodeList);
    }

    virtual ~Physics() {
        boundaries.clear();
    }

    template <typename T>
    void
    EnrollFields(std::initializer_list<const std::string> fields) {
        for (const std::string& name : fields) {
            nodeList->insertField<T>(name);
        }
    }

    template <typename T>
    void
    EnrollStateFields(std::initializer_list<const std::string> fields) {
        for (const std::string& name : fields) {
            Field<T>* field = nodeList->getField<T>(name);
            state.template addField<T>(field);
        }
    }

    virtual void
    ZeroTimeInitialize() {
        UpdateState();
        InitializeBoundaries();
    }

    virtual void
    InitializeBoundaries() {
        for (Boundary<dim>* boundary : boundaries)
            boundary->ZeroTimeInitialize();
    }

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt)  {  }

    virtual void
    PreStepInitialize() {
        state.updateLastDt(lastDt);
    };

    virtual void
    FinalizeStep(const State<dim>* finalState) {
        for (int i = 0; i < finalState->count(); ++i) {
            FieldBase* FieldToCopy = finalState->getFieldByIndex(i);
            Name fname = FieldToCopy->getName();
            if (auto* resultDouble = dynamic_cast<Field<double>*>(FieldToCopy)) {
                ScalarField* nField = nodeList->template getField<double>(fname.name());
                nField->copyValues(resultDouble);
            } else if (auto* resultVector = dynamic_cast<Field<Lin::Vector<dim>>*>(FieldToCopy)) {
                VectorField* nField = nodeList->template getField<Vector>(fname.name());
                nField->copyValues(resultVector);
            }
        } 
        FinalChecks();
    };

    virtual void
    FinalChecks() {};

    virtual void
    UpdateState() { state.updateFields(nodeList);}; 
    // you should not need to UpdateState() inside an integrator or physics package
    // since the controller calls it at the end of every step, however, if you roll
    // your own controller, this could get tricky.

    virtual void
    VerifyFields(NodeList* nodeList) {
        int numNodes = nodeList->size();
        if (nodeList->mass() == nullptr)
            nodeList->insertField<double>("mass"); // Add the mass field to the nodeList

        for (const std::string& name : 
            {"position", "velocity"}) {
            if (nodeList->getField<Vector>(name) == nullptr)
                nodeList->insertField<Vector>(name);
        }
    }

    virtual double 
    EstimateTimestep() const { return 0; }

    virtual NodeList*
    getNodeList() const { return nodeList; }

    virtual const State<dim>* 
    getState() const { return &state; }

    virtual double
    lastStep() const {return lastDt; }

    virtual void
    addBoundary(Boundary<dim>* boundary){
        boundaries.push_back(boundary);
    }

    virtual void
    ApplyBoundaries(State<dim>* bState) {
        if(boundaries.size() > 0)
            for (const auto& boundary : boundaries)
                boundary->ApplyBoundaries(bState,nodeList);
    }
    
    virtual std::string
    name() const { return "Physics"; }

    virtual std::string
    description() const { return "The base class for all physics packages"; }
};

#endif //PHYSICS_HH