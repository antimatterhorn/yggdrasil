// Copyright (C) 2026  Cody Raskin

#pragma once
#include "aleMeshBoundary.hh"

// Zero-gradient open boundary: the exterior state is just the interior state,
// unchanged. Simpler than OutflowGridBoundary's linear extrapolation (which
// needs a "next interior neighbor along this direction" concept that doesn't
// generalize cleanly to an unstructured mesh), but sufficient for an open
// boundary and effectively free given the ghostState() interface.
template <int dim>
class OutflowALEMeshBoundary : public ALEMeshBoundary<dim> {
public:
    using Vector = Lin::Vector<dim>;

    OutflowALEMeshBoundary(Mesh::ALEMesh<dim>* mesh) : ALEMeshBoundary<dim>(mesh) {}

    virtual void
    ghostState(size_t faceIndex,
               double rhoIn, const Vector& vIn, double uIn,
               double& rhoOut, Vector& vOut, double& uOut) const override {
        (void)faceIndex;
        rhoOut = rhoIn;
        vOut = vIn;
        uOut = uIn;
    }
};
