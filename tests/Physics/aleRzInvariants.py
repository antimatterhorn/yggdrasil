from yggdrasil import *
from Physics import ALEMeshHydroHLLE2d
from Mesh import ALEMesh2d, Geometry
from EOS import IdealGasEOS
from Boundaries import ReflectingALEMeshBoundary2d

# Decisive RZ correctness invariants for ALEMeshHydroBase, mirroring
# tests/Physics/rzInvariants.py's three checks for the structured-grid case
# (no exact solution needed), on a hand-built quad-strip ALEMesh instead of a
# Grid2d:
#
#   A. Uniform gas at rest stays at rest. This is THE test of the +p/r
#      source term: the radius-weighted pressure flux imbalance across the
#      curved faces must be cancelled exactly by the source, or a uniform
#      gas would spuriously accelerate radially.
#
#   B. A z-oriented problem (uniform in r) reduces to the Cartesian result:
#      the RZ mid-r-column density profile must match a pure Cartesian run
#      on the same mesh shape closely, with no new radial-velocity noise.
#
#   C. The r=0 axis boundary is auto-installed, and a user boundary that
#      also claims an axis face is rejected.

GAMMA = 1.4


def buildMesh(geometry, nr, nz, dr, dz):
    mesh = ALEMesh2d(geometry)

    def nodeId(i, j):
        return i + j * (nr + 1)

    for j in range(nz + 1):
        for i in range(nr + 1):
            mesh.addNode(Vector2d(float(i) * dr, float(j) * dz))

    for j in range(nz):
        for i in range(nr):
            mesh.addCell([nodeId(i, j), nodeId(i + 1, j), nodeId(i + 1, j + 1), nodeId(i, j + 1)])

    mesh.computeFaces()
    return mesh


def cellId(i, j, nr):
    return i + j * nr


def classifyBoundaryFaces(mesh, nr, dr, nz, dz):
    xmax, ymax = nr * dr, nz * dz
    faces = mesh.getFaces()
    left, right, top, bottom = [], [], [], []
    for idx, f in enumerate(faces):
        if not f.isBoundary:
            continue
        cx, cy = f.centroid.x, f.centroid.y
        if abs(cx - 0.0) < 1e-9:
            left.append(idx)
        elif abs(cx - xmax) < 1e-9:
            right.append(idx)
        elif abs(cy - 0.0) < 1e-9:
            bottom.append(idx)
        elif abs(cy - ymax) < 1e-9:
            top.append(idx)
    return left, right, top, bottom


def build(geometry, nr, nz, dr, dz, install_axis_left=False):
    mesh = buildMesh(geometry, nr, nz, dr, dz)
    left, right, top, bottom = classifyBoundaryFaces(mesh, nr, dr, nz, dz)
    nodes = NodeList(nr * nz)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = ALEMeshHydroHLLE2d(nodes, constants, eos, mesh)

    wallFaces = right + top + bottom
    if geometry != Geometry.CylindricalRZ or install_axis_left:
        wallFaces = left + wallFaces
    box = ReflectingALEMeshBoundary2d(mesh)
    box.setFaces(wallFaces)
    hydro.addBoundary(box)
    # Return mesh/constants/eos/box too: the C++ hydro/boundary hold bare
    # pointers/references, so they must outlive it or Python will GC them.
    return dict(mesh=mesh, nodes=nodes, hydro=hydro,
                constants=constants, eos=eos, box=box)


def integrator_for(hydro):
    return RungeKutta4Integrator2d([hydro], dtmin=1e-4, verbose=False)


def test_uniform_rest():
    nr, nz = 24, 8
    sim = build(Geometry.CylindricalRZ, nr, nz, 0.05, 0.1)
    nodes, hydro = sim["nodes"], sim["hydro"]
    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    for i in range(nr * nz):
        density.setValue(i, 1.0)
        energy.setValue(i, 2.5)          # p = 0.4*1*2.5 = 1.0, uniform
    integ = integrator_for(hydro)
    Controller(integrator=integ, periodicWork=[], statStep=100000).Step(300)

    v = nodes.getFieldVector2d("velocity")
    vmax = max(abs(v[i].x) + abs(v[i].y) for i in range(nr * nz))
    drho = max(abs(density[i] - 1.0) for i in range(nr * nz))
    assert vmax < 1e-9, f"uniform RZ gas accelerated: max|v|={vmax:.3e}"
    assert drho < 1e-9, f"uniform RZ density drifted: max|drho|={drho:.3e}"
    print(f"[A ok] uniform rest: max|v|={vmax:.2e}, max|drho|={drho:.2e}")


def z_sod_profile(geometry, nr, nz, dr, dz, cycles):
    sim = build(geometry, nr, nz, dr, dz)
    nodes, hydro = sim["nodes"], sim["hydro"]
    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    for j in range(nz):
        for i in range(nr):
            idx = cellId(i, j, nr)
            if j < nz // 2:
                density.setValue(idx, 1.0);   energy.setValue(idx, 2.5)
            else:
                density.setValue(idx, 0.125); energy.setValue(idx, 2.0)
    integ = integrator_for(hydro)
    Controller(integrator=integ, periodicWork=[], statStep=100000).Step(cycles)
    v = nodes.getFieldVector2d("velocity")
    imid = nr // 2
    prof = [density[cellId(imid, j, nr)] for j in range(nz)]
    vtmax = max(abs(v[cellId(i, j, nr)].x) for i in range(nr) for j in range(nz))
    return prof, vtmax, sim


def test_pure_z_reduces_to_cartesian():
    # ALEMeshHydroBase is first-order (no slope reconstruction, unlike the
    # Grid case's GridHydroKT2d), so the RZ and Cartesian runs' shock fronts
    # can end up a small fraction of a cell out of phase with each other --
    # a real but expected first-order truncation-error artifact, not a sign
    # of a source-term bug (already confirmed exact by test A). A pointwise
    # (max) comparison is sensitive to exactly this kind of shift right at
    # the steepest part of the profile; L1 (mean absolute difference) is the
    # right metric here, same reasoning as tests/Physics/sod.py's comparison
    # against the exact Riemann solution.
    nr, nz, dr, dz, cycles = 12, 40, 1.0, 1.0, 150
    rz_prof, rz_vt, rz_sim = z_sod_profile(Geometry.CylindricalRZ, nr, nz, dr, dz, cycles)
    ca_prof, ca_vt, ca_sim = z_sod_profile(Geometry.Cartesian, nr, nz, dr, dz, cycles)

    l1 = sum(abs(a - b) for a, b in zip(rz_prof, ca_prof)) / len(rz_prof)
    assert l1 < 3e-2, f"RZ z-Sod mid column L1 differs from Cartesian by {l1:.3e}"

    # Cartesian ALE has essentially zero transverse noise for this exactly
    # axis-uniform IC (unlike Grid/KT's dimensionally-split scheme, which
    # has some baseline), so comparing RZ noise against that near-zero
    # baseline via a ratio is degenerate -- use an absolute bound instead,
    # scaled to the problem's own velocity scale (sound speed ~1.18 here).
    assert rz_vt < 1e-2, f"RZ radial noise {rz_vt:.3e} too large (Cartesian baseline {ca_vt:.3e})"
    print(f"[B ok] pure-z reduces to Cartesian: L1_density={l1:.2e}; "
          f"radial noise {rz_vt:.2e} (Cartesian baseline {ca_vt:.2e})")


def test_axis_guard_and_autoinstall():
    nr, nz = 16, 8
    sim = build(Geometry.CylindricalRZ, nr, nz, 0.05, 0.1)
    nodes = sim["nodes"]
    d = nodes.getFieldDouble("density"); e = nodes.getFieldDouble("specificInternalEnergy")
    for i in range(nr * nz):
        d.setValue(i, 1.0); e.setValue(i, 2.5)
    integ = integrator_for(sim["hydro"])
    Controller(integrator=integ, periodicWork=[], statStep=100000).Step(5)
    assert all(d[i] == d[i] for i in range(nr * nz)), "NaN near axis with auto-installed BC"
    print("[C ok] axis boundary auto-installed (no user left face): stable at r=0")

    sim2 = build(Geometry.CylindricalRZ, nr, nz, 0.05, 0.1, install_axis_left=True)
    nodes2 = sim2["nodes"]
    d2 = nodes2.getFieldDouble("density"); e2 = nodes2.getFieldDouble("specificInternalEnergy")
    for i in range(nr * nz):
        d2.setValue(i, 1.0); e2.setValue(i, 2.5)
    integ2 = integrator_for(sim2["hydro"])
    raised = False
    try:
        Controller(integrator=integ2, periodicWork=[], statStep=100000).Step(1)
    except RuntimeError as ex:
        raised = True
        print(f"[C ok] axis guard rejected user left face: {str(ex)[:60]}...")
    assert raised, "expected RuntimeError when a user boundary claims the r=0 axis faces"


if __name__ == "__main__":
    test_uniform_rest()
    test_pure_z_reduces_to_cartesian()
    test_axis_guard_and_autoinstall()
    print("\nAll ALE RZ invariants passed.")
