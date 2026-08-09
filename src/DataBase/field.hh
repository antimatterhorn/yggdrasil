// Copyright (C) 2026  Cody Raskin

#ifndef FIELD_HH
#define FIELD_HH

#include <string>
#include <vector>
#include <memory>
#include <complex>
#include "../Type/name.hh"
#include "../Math/vectorMath.hh"

// Identifies the element type stored in a Field<T> without needing to know T,
// so generic code (e.g. restart serialization) can dispatch by runtime tag
// instead of guessing which Field<T> a FieldBase* actually points to.
enum class FieldDataType {
    Int,
    Float,
    Double,
    Complex,
    Vector1,
    Vector2,
    Vector3,
    String,
    Unknown
};

template <typename T>
struct FieldDataTypeTraits {
    static constexpr FieldDataType value = FieldDataType::Unknown;
};
template <> struct FieldDataTypeTraits<int> { static constexpr FieldDataType value = FieldDataType::Int; };
template <> struct FieldDataTypeTraits<float> { static constexpr FieldDataType value = FieldDataType::Float; };
template <> struct FieldDataTypeTraits<double> { static constexpr FieldDataType value = FieldDataType::Double; };
template <> struct FieldDataTypeTraits<std::complex<double>> { static constexpr FieldDataType value = FieldDataType::Complex; };
template <> struct FieldDataTypeTraits<std::string> { static constexpr FieldDataType value = FieldDataType::String; };
template <> struct FieldDataTypeTraits<Lin::Vector<1>> { static constexpr FieldDataType value = FieldDataType::Vector1; };
template <> struct FieldDataTypeTraits<Lin::Vector<2>> { static constexpr FieldDataType value = FieldDataType::Vector2; };
template <> struct FieldDataTypeTraits<Lin::Vector<3>> { static constexpr FieldDataType value = FieldDataType::Vector3; };

// Base class for all Field types
class FieldBase {
public:
    virtual ~FieldBase() {}  // Make the base class polymorphic

    virtual bool hasName() const = 0;
    virtual Name getName() const = 0;
    virtual std::string getNameString() const = 0;
    virtual unsigned int getSize() const = 0;
    virtual unsigned int size() const = 0;
    virtual FieldDataType type() const = 0;
};

template <typename T>
class Field : public FieldBase {
private:
    std::vector<T> values;
    Name name;

public:
    using ElementType = T;
    
    Field();
    Field(const std::string& fieldName);
    Field(const std::string& fieldName, unsigned int numElements);

    void addValue(const T& value);
    unsigned int getSize() const override;
    unsigned int size() const override;
    const std::vector<T>& getValues() const;
    const T& getValue(const unsigned int index) const;
    void setValue(const unsigned int index, T val);
    T& operator[](const unsigned int index);
    const T& operator[](const unsigned int index) const;

    void copyValues(const Field<T>& other);
    void copyValues(const Field<T>* other);

    Field<T>& operator=(const Field<T>& other);
    Field<T> operator+(const Field<T>& other) const;
    Field<T> operator-(const Field<T>& other) const;
    Field<T> operator*(const double other) const;
    Field<T>& operator+=(const Field<T>& other);
    Field<T>& operator-=(const Field<T>& other);
    Field<T>& operator*=(const double other);

    // Fused this += scalar*other, so a scaled add never builds a scaled copy first.
    Field<T>& axpy(const double scalar, const Field<T>& other);

    bool hasName() const override;
    Name getName() const override;
    std::string getNameString() const override;

    FieldDataType type() const override { return FieldDataTypeTraits<T>::value; }

    void fill(unsigned int n, T val);

    inline
    T* data() {
        return values.data();
    }

    inline
    const T* data() const {
        return values.data();
    }
};

#include "field.cc"

#endif // FIELD_HH
