// Copyright (C) 2026  Cody Raskin

#ifndef FIELD_CC
#define FIELD_CC

#include "field.hh"
#include <cstdlib>

// Below this element count an OpenMP fork/join costs more than the loop saves.
static constexpr int fieldParallelThreshold = 4096;

template <typename T>
Field<T>::Field() {}

template <typename T>
Field<T>::Field(const std::string& fieldName)
    : name(fieldName) {}

template <typename T>
Field<T>::Field(const std::string& fieldName, unsigned int numElements)
    : name(fieldName) {
    this->fill(numElements, T());
}

template <typename T>
void 
Field<T>::addValue(const T& value) {
    values.push_back(value);
}

template <typename T>
unsigned int 
Field<T>::getSize() const { 
    return values.size(); 
}

template <typename T>
unsigned int 
Field<T>::size() const { 
    return values.size(); 
}

template <typename T>
const std::vector<T>& 
Field<T>::getValues() const { 
    return values; 
}

template <typename T>
const T& 
Field<T>::getValue(const unsigned int index) const {
    if (index >= values.size()) {
        std::cerr << "you've requested item " << index << " out of " << values.size() <<
            " in " << this->getNameString() << " which doesn't exist" << std::endl;
        std::exit(1);
    }
    return values[index];
}

template <typename T>
void 
Field<T>::setValue(const unsigned int index, T val) { 
    values[index] = val; 
}

template <typename T>
T& 
Field<T>::operator[](const unsigned int index) { 
    return values[index]; 
}

template <typename T>
const T& 
Field<T>::operator[](const unsigned int index) const { 
    return values[index]; 
}

template <typename T> 
void 
Field<T>::copyValues(const Field<T>& other) {
    if (this != &other) 
        values = other.values;
}

template <typename T> 
void
Field<T>::copyValues(const Field<T>* other) {
    values = other->values;
}

template <typename T> 
Field<T>& 
Field<T>::operator=(const Field<T>& other) {
    if (this != &other) { // Avoid self-assignment
        values = other.values;
        name = other.name;
    }
    return *this;
}

template <typename T>
Field<T>
Field<T>::operator+(const Field<T>& other) const {
    Field<T> result(*this); // Create a copy of the current object
    result += other;
    return result; // Return the result
}

template <typename T>
Field<T>
Field<T>::operator-(const Field<T>& other) const {
    Field<T> result(*this); // Create a copy of the current object
    result -= other;
    return result; // Return the result
}

template <typename T>
Field<T>
Field<T>::operator*(const double other) const {
    Field<T> result(*this); // Create a copy of the current object
    result *= other;
    return result; // Return the result
}

template <typename T>
Field<T>&
Field<T>::operator+=(const Field<T>& other) {
    if (this != &other) {
        const int n = (int)values.size();
        T* a = values.data();
        const T* b = other.values.data();
        #pragma omp parallel for if(n >= fieldParallelThreshold)
        for (int i = 0; i < n; ++i) a[i] += b[i];
    }
    return *this;
}

template <typename T>
Field<T>&
Field<T>::operator-=(const Field<T>& other) {
    if (this->size() != other.size())
        throw std::invalid_argument("Field sizes do not match for subtraction");
    const int n = (int)values.size();
    T* a = values.data();
    const T* b = other.values.data();
    #pragma omp parallel for if(n >= fieldParallelThreshold)
    for (int i = 0; i < n; ++i) a[i] -= b[i];
    return *this;
}


template <typename T>
Field<T>&
Field<T>::operator*=(const double other) {
    const int n = (int)values.size();
    T* a = values.data();
    #pragma omp parallel for if(n >= fieldParallelThreshold)
    for (int i = 0; i < n; ++i) a[i] = a[i] * other;
    return *this;
}

// Lin::Vector has operator* and operator+= but no operator*=, hence this form.
template <typename T>
Field<T>&
Field<T>::axpy(const double scalar, const Field<T>& other) {
    if (this->size() != other.size())
        throw std::invalid_argument("Field sizes do not match for axpy");
    const int n = (int)values.size();
    T* a = values.data();
    const T* b = other.values.data();
    #pragma omp parallel for if(n >= fieldParallelThreshold)
    for (int i = 0; i < n; ++i) a[i] += b[i] * scalar;
    return *this;
}

template <typename T>
bool 
Field<T>::hasName() const {
    return !name.name().empty();
}

template <typename T>
Name 
Field<T>::getName() const {
    return name;
}

template <typename T>
std::string 
Field<T>::getNameString() const {
    return name.name();
}

template <typename T>
void 
Field<T>::fill(unsigned int n, T val) {
    for (unsigned int i = 0; i < n; ++i) {
        this->addValue(val);
    }
}

#endif // FIELD_CC
