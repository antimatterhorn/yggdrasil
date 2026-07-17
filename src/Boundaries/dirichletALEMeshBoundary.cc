// Copyright (C) 2026  Cody Raskin

#pragma once
#include "aleMeshBoundary.hh"

// Fixed exterior state (pinned inflow): returns the same (rho, v, u) at
// every governed face regardless of the interior state. Same role as
// DirichletGridBoundary for pinned inflow.
template <int dim>
class DirichletALEMeshBoundary : public ALEMeshBoundary<dim> {
protected:
    double rho_;
    Lin::Vector<dim> v_;
    double u_;

public:
    using Vector = Lin::Vector<dim>;

    DirichletALEMeshBoundary(Mesh::ALEMesh<dim>* mesh, double rho, Vector v, double u)
        : ALEMeshBoundary<dim>(mesh), rho_(rho), v_(v), u_(u) {}

    virtual void
    ghostState(size_t faceIndex,
               double rhoIn, const Vector& vIn, double uIn,
               double& rhoOut, Vector& vOut, double& uOut) const override {
        (void)faceIndex; (void)rhoIn; (void)vIn; (void)uIn;
        rhoOut = rho_;
        vOut = v_;
        uOut = u_;
    }
};
