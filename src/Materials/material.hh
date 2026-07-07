// Copyright (C) 2026  Cody Raskin

#pragma once

#include <Eigen/Dense>
#include <string>

enum class PlaneCondition { Stress, Strain };

template <int dim>
class FEMMaterial {
public:
    virtual ~FEMMaterial() = default;

    virtual Eigen::MatrixXd materialMatrix() const = 0;

    virtual std::string name() const = 0;
};
