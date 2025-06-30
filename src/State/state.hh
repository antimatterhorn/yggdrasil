// Copyright (C) 2025  Cody Raskin

#ifndef STATE_HH
#define STATE_HH

#include "../DataBase/nodeList.hh"
#include "../Type/name.hh"
#include "../Math/vectorMath.hh"

template <int dim>
class State {
private:   
    std::vector<std::shared_ptr<FieldBase>> fields;
    int numNodes;
    double lastDt;
public:
    using Complex       = std::complex<double>;
    using Vector        = Lin::Vector<dim>;
    using ComplexField  = Field<Complex>;
    using VectorField   = Field<Vector>;
    using ScalarField   = Field<double>;

    State(int numNodes) : 
        numNodes(numNodes) { };

    ~State() {};

    template <typename T>
    Field<T>* getFieldByName(const Name& name) const {
        for (const auto& fieldPtr : fields) {
            if (fieldPtr->hasName() && fieldPtr->getNameString() == name.name()) {
                Field<T>* castedField = dynamic_cast<Field<T>*>(fieldPtr.get());
                if (castedField != nullptr) {
                    return castedField; // Return the field if found and correctly casted
                }
            }
        }
        std::cout << "No field in this State with name " << name.name() << std::endl;
        return nullptr; // Return nullptr if no matching field is found
    }

    template <typename T>
    void addField(const Field<T>* fieldPtr) {
        Name name = fieldPtr->getName();
        std::shared_ptr<Field<T>> newField = std::make_shared<Field<T>>(name.name(), this->size());
        newField->copyValues(fieldPtr);
        fields.push_back(newField);
    }

    template <typename T>
    Field<T>* 
    getField(const std::string& name) const {
        return getFieldByName<T>(Name(name));
    }
    
    template <typename T>
    void 
    insertField(const std::string& name) {
        auto newField = std::make_shared<Field<T>>(name, this->size());
        fields.push_back(newField); // Use make_shared for field creation
    }

    FieldBase* 
    getFieldByIndex(int index) const {
        if (index < 0 || index >= fields.size()) {
            return nullptr; // Return nullptr if index is out of range
        }
        return fields[index].get(); // Return the raw pointer to the field at index
    }

    void
    updateFields(NodeList* nodeList) {
        for (const auto& fieldPtr : fields) {
            if (auto* doubleField = dynamic_cast<ScalarField*>(fieldPtr.get())) {
                auto* otherDoubleField = dynamic_cast<const ScalarField*>(nodeList->getField<double>(doubleField->getNameString()));
                if (otherDoubleField) {
                    doubleField->copyValues(otherDoubleField);
                }
            } else if (auto* vectorField = dynamic_cast<VectorField*>(fieldPtr.get())) {
                auto* otherVectorField = dynamic_cast<const VectorField*>(nodeList->getField<Vector>(vectorField->getNameString()));
                if (otherVectorField) {
                    vectorField->copyValues(otherVectorField);
                }
            } else if (auto* complexField = dynamic_cast<ComplexField*>(fieldPtr.get())) {
                auto* otherComplexField = dynamic_cast<const ComplexField*>(nodeList->getField<Complex>(complexField->getNameString()));
                if (otherComplexField) {
                    complexField->copyValues(otherComplexField);
                }
            }
        }
    }

    int
    size() const { return numNodes; }

    int
    count() const { return fields.size(); }

    void
    addState(const State* other) {
        if (this->count() != other->count() || this->size() != other->size()) {
            std::cout << this->count() << other->count() << this->size() << other->size() << std::endl;
            throw std::invalid_argument("Incompatible State objects for addition");
        }
        for (const auto& fieldPtr : fields) {
            if (auto* doubleField = dynamic_cast<ScalarField*>(fieldPtr.get())) {
                auto* otherDoubleField = dynamic_cast<const ScalarField*>(other->getField<double>(doubleField->getNameString()));
                if (otherDoubleField) {
                    *doubleField += *otherDoubleField;
                }
            } else if (auto* vectorField = dynamic_cast<VectorField*>(fieldPtr.get())) {
                auto* otherVectorField = dynamic_cast<const VectorField*>(other->getField<Vector>(vectorField->getNameString()));
                if (otherVectorField) {
                    *vectorField += *otherVectorField;
                }
            } else if (auto* complexField = dynamic_cast<ComplexField*>(fieldPtr.get())) {
                auto* otherComplexField = dynamic_cast<const ComplexField*>(other->getField<Complex>(complexField->getNameString()));
                if (otherComplexField) {
                    *complexField += *otherComplexField;
                }
            }
        }
    }

    State& 
    operator+=(const State& other) {
        if (this != &other) {
            if (this->count() != other.count() || this->size() != other.size()) {
                std::cout << this->count() << other.count() << this->size() << other.size() << std::endl;
                throw std::invalid_argument("Incompatible State objects for addition");
            }
            for (const auto& fieldPtr : fields) {
                if (auto* doubleField = dynamic_cast<ScalarField*>(fieldPtr.get())) {
                    auto* otherDoubleField = dynamic_cast<const ScalarField*>(other.getField<double>(doubleField->getNameString()));
                    if (otherDoubleField) {
                        *doubleField += *otherDoubleField;
                    }
                } else if (auto* vectorField = dynamic_cast<VectorField*>(fieldPtr.get())) {
                    auto* otherVectorField = dynamic_cast<const VectorField*>(other.getField<Vector>(vectorField->getNameString()));
                    if (otherVectorField) {
                        *vectorField += *otherVectorField;
                    }
                } else if (auto* complexField = dynamic_cast<ComplexField*>(fieldPtr.get())) {
                    auto* otherComplexField = dynamic_cast<const ComplexField*>(other.getField<Complex>(complexField->getNameString()));
                    if (otherComplexField) {
                        *complexField += *otherComplexField;
                    }
                }
            }
        }
        return *this;
    }

    // Implement *= operator for State class
    State& 
    operator*=(const double other) {
        for (const auto& fieldPtr : fields) {
            if (auto* doubleField = dynamic_cast<ScalarField*>(fieldPtr.get())) {
                *doubleField *= other;
            } else if (auto* vectorField = dynamic_cast<VectorField*>(fieldPtr.get())) {
                *vectorField *= other;
            } else if (auto* complexField = dynamic_cast<ComplexField*>(fieldPtr.get()))  {
                *complexField *= other;
            }
        }
        return *this;
    }

    State<dim> operator*(const double& scalar) const {
        State<dim> newState(numNodes); // Create a new State with the same number of nodes
        newState.clone(this); // Clone the fields from the current state to the new state
        
        for (int i = 0; i < newState.count(); ++i) {
            FieldBase* fieldPtr = newState.getFieldByIndex(i);
            if (auto* doubleField = dynamic_cast<ScalarField*>(fieldPtr)) {
                *doubleField *= scalar;
            } else if (auto* vectorField = dynamic_cast<VectorField*>(fieldPtr)) {
                *vectorField *= scalar;
            } else if (auto* complexField = dynamic_cast<ComplexField*>(fieldPtr))  {
                *complexField *= scalar;
            }
        }
        return newState;
    }

    // Overload the = assignment operator to copy the fields from another state
    State& 
    operator=(const State& rhs) {
        fields = rhs.fields;
        numNodes = rhs.numNodes;
        return *this;
    }

    State 
    operator+(const State& other) const {
        if (this->count() != other.count() || this->size() != other.size()) {
            std::cout << this->count() << other.count() << this->size() << other.size() << std::endl;
            throw std::invalid_argument("Incompatible State objects for addition");
        }

        for (int i = 0; i < this->count(); ++i) {
            const auto& fieldPtr = this->getFieldByIndex(i);
            const auto& otherFieldPtr = other.getFieldByIndex(i);
            
            if (fieldPtr->getNameString() != otherFieldPtr->getNameString()) {
                throw std::invalid_argument("Incompatible State objects for addition: mismatched field names");
            }
            if (typeid(*fieldPtr) != typeid(*otherFieldPtr)) {
                throw std::invalid_argument("Incompatible State objects for addition: mismatched field types");
            }
        }

        State<dim> result(this->size()); // Create a new State object with the same size

        for (int i = 0; i < this->count(); ++i) {
            FieldBase* field = this->getFieldByIndex(i); // Get the field at index i

            if (auto* doubleField = dynamic_cast<ScalarField*>(field)) {
                auto* otherDoubleField = dynamic_cast<const ScalarField*>(other.getFieldByIndex(i));
                if (otherDoubleField) {
                    *doubleField += *otherDoubleField; // Perform addition for double fields
                }
            } else if (auto* vectorField = dynamic_cast<VectorField*>(field)) {
                auto* otherVectorField = dynamic_cast<const VectorField*>(other.getFieldByIndex(i));
                if (otherVectorField) {
                    *vectorField += *otherVectorField; // Perform addition for vector fields
                }
            } else if (auto* complexField = dynamic_cast<ComplexField*>(field)) {
                auto* otherComplexField = dynamic_cast<const ComplexField*>(other.getFieldByIndex(i));
                if (otherComplexField) {
                    *complexField += *otherComplexField; // Perform addition for vector fields
                }
            }
        }

        return result;
    }

    State<dim> operator-(const State& other) const {
        if (this->count() != other.count() || this->size() != other.size()) {
            throw std::invalid_argument("Incompatible State objects for subtraction");
        }

        State<dim> result(this->size());
        result.clone(this);

        for (int i = 0; i < result.count(); ++i) {
            FieldBase* resultField = result.getFieldByIndex(i);
            FieldBase* otherField = other.getFieldByIndex(i);

            if (auto* resultDouble = dynamic_cast<ScalarField*>(resultField)) {
                auto* otherDouble = dynamic_cast<const ScalarField*>(otherField);
                if (otherDouble) *resultDouble -= *otherDouble;
            } else if (auto* resultVector = dynamic_cast<VectorField*>(resultField)) {
                auto* otherVector = dynamic_cast<const VectorField*>(otherField);
                if (otherVector) *resultVector -= *otherVector;
            } else if (auto* resultComplex = dynamic_cast<ComplexField*>(resultField)) {
                auto* otherComplex = dynamic_cast<const ComplexField*>(otherField);
                if (otherComplex) *resultComplex -= *otherComplex;
            }
        }

        return result;
    }


    void
    clone(const State* other) {
        fields.clear(); 

        for (int i = 0; i < other->count(); ++i) {
            FieldBase* field = other->getFieldByIndex(i);
            if (field->hasName()) {
                std::string fieldName = field->getNameString();
                if (dynamic_cast<ScalarField*>(field) != nullptr) {
                    insertField<double>(fieldName);
                    ScalarField* thisField = getField<double>(fieldName);
                    thisField->copyValues(dynamic_cast<ScalarField*>(field));
                } 
                else if (dynamic_cast<VectorField*>(field) != nullptr) {
                    insertField<Vector>(fieldName);
                    VectorField* thisField = getField<Vector>(fieldName);
                    thisField->copyValues(dynamic_cast<VectorField*>(field));
                }
                else if (dynamic_cast<ComplexField*>(field) != nullptr) {
                    insertField<Complex>(fieldName);
                    ComplexField* thisField = getField<Complex>(fieldName);
                    thisField->copyValues(dynamic_cast<ComplexField*>(field));
                }
            }
        }

    }

    State<dim> deepCopy() const {
        State<dim> out(this->size());
        out.clone(this);
        return out;
    }

    void
    ghost(const State* other) {
        fields.clear(); 
        
        for (int i = 0; i < other->count(); ++i) {
            FieldBase* field = other->getFieldByIndex(i);
            if (field->hasName()) {
                std::string fieldName = field->getNameString();
                if (dynamic_cast<ScalarField*>(field) != nullptr) {
                    insertField<double>(fieldName);
                } 
                else if (dynamic_cast<VectorField*>(field) != nullptr) {
                    insertField<Vector>(fieldName);
                }
                else if (dynamic_cast<ComplexField*>(field) != nullptr) {
                    insertField<Complex>(fieldName);
                }
            }
        }

    }

    double L2Norm() const {
        double sum = 0.0;

        for (int i = 0; i < this->count(); ++i) {
            FieldBase* field = this->getFieldByIndex(i);

            if (auto* f = dynamic_cast<ScalarField*>(field)) {
                for (int j = 0; j < f->size(); ++j)
                    sum += (*f)[j] * (*f)[j];
            } else if (auto* f = dynamic_cast<VectorField*>(field)) {
                for (int j = 0; j < f->size(); ++j)
                    sum += (*f)[j].mag2();
            } else if (auto* f = dynamic_cast<ComplexField*>(field)) {
                for (int j = 0; j < f->size(); ++j)
                    sum += std::norm((*f)[j]);
            }
        }

        return std::sqrt(sum);
    }

    void swap(State& other) {
        std::swap(this->fields, other.fields);
        std::swap(this->numNodes, other.numNodes);
        std::swap(this->lastDt, other.lastDt);
    }


    void
    updateLastDt(const double dt) { lastDt = dt;}

    double
    getLastDt() { return lastDt;}
};


#endif //STATE_HH