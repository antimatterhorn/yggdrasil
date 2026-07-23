// Copyright (C) 2026  Cody Raskin

#pragma once

#include "boundary.hh"
#include "../Mesh/alemesh.hh"

// Base for ALEMesh boundary conditions. Unlike GridBoundary (which mutates
// ghost-cell field values via ApplyScalar/ApplyVector/ApplyComplex), an
// ALEMesh boundary condition answers a different question: "given the real
// interior state at this boundary face, what exterior/ghost state should the
// Riemann solver pair it against?" ALEMeshHydroBase calls ghostState()
// directly for each boundary face each step, rather than going through the
// generic Boundary::ApplyBoundaries mechanism (same reasoning that already
// motivated Boundary::getObstacleIds() as a purpose-built extension point).
template <int dim>
class ALEMeshBoundary : public Boundary<dim> {
protected:
    Mesh::ALEMesh<dim>* mesh;
    std::vector<size_t> faceIds;

public:
    using Vector = Lin::Vector<dim>;

    ALEMeshBoundary(Mesh::ALEMesh<dim>* mesh) : mesh(mesh) {}
    virtual ~ALEMeshBoundary() {}

    // Explicitly register which boundary faces this condition governs.
    virtual void setFaces(const std::vector<size_t>& ids) { faceIds = ids; }

    // Convenience: register every boundary face currently on the mesh
    // (Face::isBoundary()). Call after mesh->computeFaces().
    virtual void setAllBoundaryFaces() {
        faceIds.clear();
        const auto& faces = mesh->getFaces();
        for (size_t i = 0; i < faces.size(); ++i)
            if (faces[i].isBoundary()) faceIds.push_back(i);
    }

    const std::vector<size_t>& getFaceIds() const { return faceIds; }

    // EOS-agnostic: only primitive density/velocity/specific-internal-energy
    // in and out. ALEMeshHydroBase evaluates the EOS itself on whatever
    // (rhoOut, uOut) come back to get pressure/soundSpeed, the same way
    // GridBoundary never needs to know about EquationOfState either.
    virtual void ghostState(size_t faceIndex,
                             double rhoIn, const Vector& vIn, double uIn,
                             double& rhoOut, Vector& vOut, double& uOut) const = 0;
};
