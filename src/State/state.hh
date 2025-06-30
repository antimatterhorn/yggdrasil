// Copyright (C) 2025  Cody Raskin

#ifndef STATE_HH
#define STATE_HH

#include "../DataBase/nodeList.hh"
#include "../Type/name.hh"
#include "../Math/vectorMath.hh"

#include <complex>
#include <memory>
#include <vector>
#include <iostream>
#include <stdexcept>


template <int dim>
class State {
private:
    std::vector<std::shared_ptr<FieldBase>> fields;
    int numNodes;
    double lastDt;

    using Complex       = std::complex<double>;
    using Vector        = Lin::Vector<dim>;
    using ComplexField  = Field<Complex>;
    using VectorField   = Field<Vector>;
    using ScalarField   = Field<double>;

    template <typename Func>
    void type_dispatch(FieldBase* f, Func&& func) const {
        if (auto* fd = dynamic_cast<ScalarField*>(f)) func(fd);
        else if (auto* fv = dynamic_cast<VectorField*>(f)) func(fv);
        else if (auto* fc = dynamic_cast<ComplexField*>(f)) func(fc);
    }

    template <typename Func>
    void forEachTypedField(Func&& func) {
        for (auto& fieldPtr : fields) {
            type_dispatch(fieldPtr.get(), func);
        }
    }

    template <typename Func>
    void forEachTypedField(Func&& func) const {
        for (const auto& fieldPtr : fields) {
            type_dispatch(fieldPtr.get(), func);
        }
    }

    template <typename Func>
    void forEachMatchingTypedField(const State<dim>& other, Func&& func) {
        for (int i = 0; i < this->count(); ++i) {
            auto* a = this->getFieldByIndex(i);
            auto* b = other.getFieldByIndex(i);
            type_dispatch(a, [&](auto* fa) {
                using T = typename std::decay_t<decltype(*fa)>::ElementType;
                if (auto* fb = dynamic_cast<Field<T>*>(b)) func(fa, fb);
            });
        }
    }

public:
    State(int numNodes) : numNodes(numNodes) {}
    ~State() {}

    template <typename T>
    Field<T>* getFieldByName(const Name& name) const {
        for (const auto& fieldPtr : fields) {
            if (fieldPtr->hasName() && fieldPtr->getNameString() == name.name()) {
                return dynamic_cast<Field<T>*>(fieldPtr.get());
            }
        }
        std::cout << "No field in this State with name " << name.name() << std::endl;
        return nullptr;
    }

    template <typename T>
    Field<T>* getField(const std::string& name) const {
        return getFieldByName<T>(Name(name));
    }

    template <typename T>
    void addField(const Field<T>* fieldPtr) {
        Name name = fieldPtr->getName();
        auto newField = std::make_shared<Field<T>>(name.name(), this->size());
        newField->copyValues(fieldPtr);
        fields.push_back(newField);
    }

    template <typename T>
    void insertField(const std::string& name) {
        fields.push_back(std::make_shared<Field<T>>(name, this->size()));
    }

    FieldBase* getFieldByIndex(int index) const {
        if (index < 0 || index >= fields.size()) return nullptr;
        return fields[index].get();
    }

    void updateFields(NodeList* nodeList) {
        forEachTypedField([&](auto* f) {
            using T = typename std::decay_t<decltype(*f)>::ElementType;
            auto* source = nodeList->getField<T>(f->getNameString());
            if (source) f->copyValues(source);
        });
    }

    int size() const { return numNodes; }
    int count() const { return fields.size(); }

    void addState(const State* other) { *this += *other; }

    State& operator+=(const State& other) {
        if (this != &other) {
            if (this->count() != other.count() || this->size() != other.size()) {
                throw std::invalid_argument("Incompatible State objects for addition");
            }
            forEachMatchingTypedField(other, [](auto* a, const auto* b) { *a += *b; });
        }
        return *this;
    }

    State& operator*=(const double other) {
        forEachTypedField([&](auto* f) { *f *= other; });
        return *this;
    }

    State operator*(const double& scalar) const {
        State<dim> newState(numNodes);
        newState.clone(this);
        newState.forEachTypedField([&](auto* f) { *f *= scalar; });
        return newState;
    }

    State& operator=(const State& rhs) {
        fields = rhs.fields;
        numNodes = rhs.numNodes;
        return *this;
    }

    State operator+(const State& other) const {
        State<dim> result(this->size());
        result.clone(this);
        result += other;
        return result;
    }

    State<dim> operator-(const State<dim>& other) const {
        State<dim> result(this->size());
        result.clone(this);
        result.forEachMatchingTypedField(other, [](auto* a, const auto* b) { *a -= *b; });
        return result;
    }

    void clone(const State* other) {
        fields.clear();
        for (int i = 0; i < other->count(); ++i) {
            FieldBase* field = other->getFieldByIndex(i);
            if (!field->hasName()) continue;

            std::string name = field->getNameString();
            if (dynamic_cast<ScalarField*>(field)) {
                insertField<double>(name);
                getField<double>(name)->copyValues(dynamic_cast<ScalarField*>(field));
            } else if (dynamic_cast<VectorField*>(field)) {
                insertField<Vector>(name);
                getField<Vector>(name)->copyValues(dynamic_cast<VectorField*>(field));
            } else if (dynamic_cast<ComplexField*>(field)) {
                insertField<Complex>(name);
                getField<Complex>(name)->copyValues(dynamic_cast<ComplexField*>(field));
            }
        }
    }

    State<dim> deepCopy() const {
        State<dim> out(this->size());
        out.clone(this);
        return out;
    }

    void ghost(const State* other) {
        fields.clear();
        for (int i = 0; i < other->count(); ++i) {
            FieldBase* field = other->getFieldByIndex(i);
            if (!field->hasName()) continue;
            std::string name = field->getNameString();

            if (dynamic_cast<ScalarField*>(field)) insertField<double>(name);
            else if (dynamic_cast<VectorField*>(field)) insertField<Vector>(name);
            else if (dynamic_cast<ComplexField*>(field)) insertField<Complex>(name);
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

    void updateLastDt(const double dt) { lastDt = dt; }
    double getLastDt() { return lastDt; }
};

#endif // STATE_HH
