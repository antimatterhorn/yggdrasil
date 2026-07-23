// Copyright (C) 2026  Cody Raskin

#pragma once
#include "aleMeshHydroBase.hh"
#include "HLL.cc"

// First-order HLLC on ALEMesh: the exact contact-preserving Riemann solve,
// with no slope reconstruction (that needs face-adjacency-chain walking to
// find a "next" cell along a direction, which doesn't exist on an
// unstructured mesh yet -- see aleMeshHydroBase.hh and the ALE hydro plan).
template<int dim>
class ALEMeshHydroHLLC : public ALEMeshHydroBase<dim> {
public:
    using typename ALEMeshHydroBase<dim>::Vector;
    using typename ALEMeshHydroBase<dim>::VectorField;
    using typename ALEMeshHydroBase<dim>::ScalarField;

    ALEMeshHydroHLLC(NodeList* nodeList, PhysicalConstants& constants,
                      EquationOfState* eos, Mesh::ALEMesh<dim>* mesh,
                      NodeList* nodeVelocities = nullptr)
        : ALEMeshHydroBase<dim>(nodeList, constants, eos, mesh, nodeVelocities) {}

    virtual std::string name() const override { return "ALEMeshHydroHLLC"; }

    virtual HLLFlux<dim>
    computeFlux(double rhoL, const Vector& vL, double uL, double pL, double csL,
                double rhoR, const Vector& vR, double uR, double pR, double csR,
                const Vector& normal, double faceVelocity) const override {
        return computeHLLCFluxFromStatesALE<dim>(rhoL, vL, uL, pL, csL,
                                                  rhoR, vR, uR, pR, csR,
                                                  normal, faceVelocity);
    }
};
