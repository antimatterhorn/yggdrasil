// Copyright (C) 2026  Cody Raskin

#pragma once

namespace Mesh {
    // CylindricalRZ makes a 2D mesh axisymmetric (axis 0 = r, axis 1 = z): cell
    // volumes and face areas are weighted by radius. Shared by Grid and
    // ALEMesh -- a runtime property of the mesh's metric, not tied to either
    // mesh's own topology representation.
    enum class Geometry { Cartesian, CylindricalRZ };
}
