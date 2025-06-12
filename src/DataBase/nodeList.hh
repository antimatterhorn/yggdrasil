// Copyright (C) 2025  Cody Raskin

#ifndef NODELIST_HH
#define NODELIST_HH

#include <vector>
#include <string>
#include <memory> // Include for std::shared_ptr
#include "field.hh"
#include "../Type/name.hh"
#include "../Math/vectorMath.hh"
#include "../Mesh/element.hh"

class NodeList {
private:
    std::vector<std::shared_ptr<FieldBase>> _fields;
    Field<int> _ids;
    std::vector<std::shared_ptr<FieldBase>> _extraFields;

public:
    NodeList();
    NodeList(int numNodes);
    ~NodeList();

    void addField(const std::shared_ptr<FieldBase>& fieldPtr);

    template <typename T>
    void insertField(const std::string& name);

    void addNodeList(const NodeList& other);

    NodeList operator+(const NodeList& other) const;
    NodeList& operator+=(const NodeList& other);

    template <typename T>
    Field<T>* getFieldByName(const Name& name) const;

    template <typename T>
    Field<T>* getField(const std::string& name) const;

    Field<double>* mass() const;

    template <int dim>
    Field<Lin::Vector<dim>>* velocity() const;

    template <int dim>
    Field<Lin::Vector<dim>>* position() const;

    size_t getFieldCount() const;
    size_t getNumNodes() const;
    std::vector<std::string> fieldNames() const;
    Field<int>& nodes();
    unsigned int size() const;

    template <int dim>
    void updatePositions(const std::vector<std::array<double, dim>>& py_positions) {
        Field<Lin::Vector<dim>>& posField = *this->position<dim>();

        if (py_positions.size() != posField.size()) {
            throw std::runtime_error("updatePosition: size mismatch between input and field");
        }

        for (size_t i = 0; i < py_positions.size(); ++i) {
            Lin::Vector<dim> v;
            for (int d = 0; d < dim; ++d) {
                v[d] = py_positions[i][d];
            }
            posField.setValue(i, v);
        }
    }
};

#include "nodeList.cc" // Include the template implementation

#endif // NODELIST_HH
