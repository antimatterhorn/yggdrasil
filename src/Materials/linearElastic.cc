// Copyright (C) 2026  Cody Raskin

#include "material.hh"
#include <stdexcept>

template <int dim>
class IsotropicLinearElastic : public FEMMaterial<dim> {
    double E, nu;
    PlaneCondition condition;

public:
    IsotropicLinearElastic(double E, double nu,
                           PlaneCondition condition = PlaneCondition::Stress)
        : E(E), nu(nu), condition(condition) {}

    Eigen::MatrixXd materialMatrix() const override {
        if constexpr (dim == 2) {
            Eigen::Matrix3d D;
            if (condition == PlaneCondition::Stress) {
                double c = E / (1.0 - nu * nu);
                D << c,      c * nu,             0.0,
                     c * nu, c,                  0.0,
                     0.0,    0.0,    c * (1.0 - nu) / 2.0;
            } else {
                double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                D << c * (1.0 - nu), c * nu,         0.0,
                     c * nu,         c * (1.0 - nu), 0.0,
                     0.0,            0.0,             c * (1.0 - 2.0 * nu) / 2.0;
            }
            return D;
        } else if constexpr (dim == 3) {
            double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
            double mu     = E / (2.0 * (1.0 + nu));

            // Voigt notation: [σxx, σyy, σzz, σyz, σxz, σxy]
            Eigen::Matrix<double, 6, 6> D = Eigen::Matrix<double, 6, 6>::Zero();
            D(0,0) = lambda + 2*mu; D(0,1) = lambda;       D(0,2) = lambda;
            D(1,0) = lambda;        D(1,1) = lambda + 2*mu; D(1,2) = lambda;
            D(2,0) = lambda;        D(2,1) = lambda;        D(2,2) = lambda + 2*mu;
            D(3,3) = mu;
            D(4,4) = mu;
            D(5,5) = mu;
            return D;
        } else {
            throw std::invalid_argument("IsotropicLinearElastic: unsupported dimension");
        }
    }

    double youngsModulus() const { return E; }
    double poissonsRatio() const { return nu; }
    PlaneCondition planeCondition() const { return condition; }

    std::string name() const override { return "IsotropicLinearElastic"; }
};
