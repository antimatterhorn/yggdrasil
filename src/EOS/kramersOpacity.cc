// Copyright (C) 2026  Cody Raskin

#include "opacityModel.hh"
#include <cmath>

// Composite Rosseland-mean opacity: Kramers bound-free/free-free absorption
// (kappa ~ rho * T^-3.5) plus a constant electron-scattering floor. Standard
// textbook approximation (e.g. Kippenhahn & Weigert, "Stellar Structure and
// Evolution", ch. 17).
class KramersOpacity : public OpacityModel {
    double kappa0;   // Kramers coefficient
    double kappaES;  // electron-scattering opacity floor
public:
    KramersOpacity(double kappa0, double kappaES, PhysicalConstants& constants)
        : OpacityModel(constants), kappa0(kappa0), kappaES(kappaES) {}

    double kappa(double rho, double T) const {
        return kappa0 * rho * std::pow(T, -3.5) + kappaES;
    }

    virtual void setOpacity(Field<double>* opacity, Field<double>* density, Field<double>* temperature) const override {
        for (int i = 0; i < opacity->size(); ++i)
            opacity->setValue(i, kappa(density->getValue(i), temperature->getValue(i)));
    }
    virtual void setConductivity(Field<double>* conductivity, Field<double>* density, Field<double>* temperature) const override {
        double sig = this->constants.stefanBoltzmannConstant();
        for (int i = 0; i < conductivity->size(); ++i) {
            double T = temperature->getValue(i);
            double rho = density->getValue(i);
            conductivity->setValue(i, 16.0*sig*T*T*T/(3.0*kappa(rho,T)*rho));
        }
    }
    virtual void setOpacity(double* opacity, double* density, double* temperature) const override {
        (*opacity) = kappa(*density, *temperature);
    }
    virtual void setConductivity(double* conductivity, double* density, double* temperature) const override {
        double sig = this->constants.stefanBoltzmannConstant();
        double T = *temperature;
        double rho = *density;
        (*conductivity) = 16.0*sig*T*T*T/(3.0*kappa(rho,T)*rho);
    }

    virtual std::string name() const override { return "KramersOpacity"; }
};
