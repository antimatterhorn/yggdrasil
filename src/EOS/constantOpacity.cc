#include "opacityModel.hh"

class ConstantOpacity : public OpacityModel {
    double k0;
public:
    ConstantOpacity(double k0, PhysicalConstants& constants) : OpacityModel(constants), k0(k0)  {}

    virtual void setOpacity(Field<double>* opacity, Field<double>* density, Field<double>* temperature) const override {
        for (int i = 0 ; i < opacity->size(); ++i){
            opacity->setValue(i, k0);
        }
    }
    virtual void setConductivity(Field<double>* conductivity, Field<double>* density, Field<double>* temperature) const override {
        for (int i = 0 ; i < conductivity->size(); ++i){
            conductivity->setValue(i, 1.0/k0);
        }
    }
    virtual void setOpacity(double* opacity, double* density, double* temperature) const override {
        (*opacity) = k0;
    }
    virtual void setConductivity(double* conductivity, double* density, double* temperature) const override {
        (*conductivity) = 1.0/k0;
    }

    virtual std::string 
    name() const override { return "ConstantOpacity"; }

};