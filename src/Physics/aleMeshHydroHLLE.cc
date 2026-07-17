// Copyright (C) 2026  Cody Raskin

#pragma once
#include "aleMeshHydroBase.hh"
#include "HLL.cc"

// Two-wave HLLE on ALEMesh -- more diffusive than ALEMeshHydroHLLC (no exact
// contact resolution), but without HLLC's star-region solve, so it has no
// equivalent of the sonic-point entropy-violation weakness that limits HLLC
// to a short stable window on a strong discontinuity. First-order (no
// reconstruction), same as ALEMeshHydroHLLC -- see aleMeshHydroBase.hh.
template<int dim>
class ALEMeshHydroHLLE : public ALEMeshHydroBase<dim> {
public:
    using typename ALEMeshHydroBase<dim>::Vector;
    using typename ALEMeshHydroBase<dim>::VectorField;
    using typename ALEMeshHydroBase<dim>::ScalarField;

    ALEMeshHydroHLLE(NodeList* nodeList, PhysicalConstants& constants,
                      EquationOfState* eos, Mesh::ALEMesh<dim>* mesh)
        : ALEMeshHydroBase<dim>(nodeList, constants, eos, mesh) {}

    virtual std::string name() const override { return "ALEMeshHydroHLLE"; }

    virtual HLLFlux<dim>
    computeFlux(double rhoL, const Vector& vL, double uL, double pL, double csL,
                double rhoR, const Vector& vR, double uR, double pR, double csR,
                const Vector& normal) const override {
        return computeHLLEFluxFromStates<dim>(rhoL, vL, uL, pL, csL,
                                              rhoR, vR, uR, pR, csR,
                                              normal);
    }
};
