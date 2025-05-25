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
        double sig = this->constants.stefanBoltzmannConstant();
        for (int i = 0 ; i < conductivity->size(); ++i){
            double T = temperature->getValue(i);
            double rho = density->getValue(i);
            conductivity->setValue(i, 16.0*sig*T*T*T/(3.0*k0*rho));
        }
    }
    virtual void setOpacity(double* opacity, double* density, double* temperature) const override {
        (*opacity) = k0;
    }
    virtual void setConductivity(double* conductivity, double* density, double* temperature) const override {
        double sig = this->constants.stefanBoltzmannConstant();
        double T = *temperature;
        double rho = *density;
        (*conductivity) = 16.0*sig*T*T*T/(3.0*k0*rho);
    }

    virtual std::string 
    name() const override { return "ConstantOpacity"; }

};