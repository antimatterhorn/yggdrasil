// Copyright (C) 2026  Cody Raskin

#pragma once
#include "aleMeshBoundary.hh"

// Solid-wall boundary: mirror the normal velocity component across the face,
// zero-gradient density/energy. Same physical convention as
// ReflectingGridBoundary's face-mode vector handling.
template <int dim>
class ReflectingALEMeshBoundary : public ALEMeshBoundary<dim> {
public:
    using Vector = Lin::Vector<dim>;

    ReflectingALEMeshBoundary(Mesh::ALEMesh<dim>* mesh) : ALEMeshBoundary<dim>(mesh) {}

    virtual void
    ghostState(size_t faceIndex,
               double rhoIn, const Vector& vIn, double uIn,
               double& rhoOut, Vector& vOut, double& uOut) const override {
        const auto& face = this->mesh->getFaces()[faceIndex];
        double vn = vIn.dot(face.normal);
        vOut = vIn - face.normal * (2.0 * vn);
        rhoOut = rhoIn;
        uOut = uIn;
    }
};
