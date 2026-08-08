from yggdrasil import *
from Physics import GridHydroKT2d
from Mesh import Grid2d, Geometry
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d

# Decisive RZ correctness invariants (no exact-solution needed):
#
#   A. Uniform gas at rest stays at rest. This is THE test of the +p/r source
#      term: the radius-weighted pressure flux imbalance across the curved
#      r-faces must be cancelled exactly by the source, or a uniform gas would
#      spuriously accelerate radially.
#
#   B. A z-oriented problem (uniform in r) reduces to the Cartesian 1D-in-z
#      result: the RZ mid-column density profile must match a pure Cartesian run
#      closely, and RZ must not introduce radial velocity noise beyond what the
#      base (dimensionally split) scheme already generates transversely at a
#      grid-aligned shock in Cartesian. (A planar shock is NOT kept perfectly 1D
#      by KT/HLL in any geometry — the Cartesian solver produces the same small
#      transverse velocity — so this is a "no new RZ instability" check, not a
#      round-off check.)
#
#   C. The r=0 axis boundary is auto-installed, and a user boundary that also
#      claims the 'left' face is rejected.

GAMMA = 1.4


def build(geometry, nr, nz, dr, dz, install_axis_left=False):
    grid = Grid2d(nr, nz, dr, dz, geometry) if geometry == Geometry.CylindricalRZ \
        else Grid2d(nr, nz, dr, dz)
    nodes = NodeList(nr * nz)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroKT2d(nodes, constants, eos, grid)
    faces = ["right", "top", "bottom"]
    if geometry != Geometry.CylindricalRZ or install_axis_left:
        faces = ["left"] + faces
    box = ReflectingGridBoundary2d(grid=grid)
    box.setFaces(faces)
    hydro.addBoundary(box)
    # Return constants/eos/box too: the C++ hydro holds a bare pointer to the
    # EOS and a reference to the constants, so they must outlive it or Python
    # will GC them and leave the solver dangling.
    return dict(grid=grid, nodes=nodes, hydro=hydro,
                constants=constants, eos=eos, box=box)


def integrator_for(hydro):
    return RungeKutta4Integrator2d([hydro], dtmin=1e-4, verbose=False)


def test_uniform_rest():
    nr, nz = 24, 8
    sim = build(Geometry.CylindricalRZ, nr, nz, 0.05, 0.1)
    grid, nodes, hydro = sim["grid"], sim["nodes"], sim["hydro"]
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
    grid, nodes, hydro = sim["grid"], sim["nodes"], sim["hydro"]
    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    for j in range(nz):
        for i in range(nr):
            idx = grid.index(i, j, 0)
            if j < nz // 2:
                density.setValue(idx, 1.0);   energy.setValue(idx, 2.5)
            else:
                density.setValue(idx, 0.125); energy.setValue(idx, 2.0)
    integ = integrator_for(hydro)
    Controller(integrator=integ, periodicWork=[], statStep=100000).Step(cycles)
    # mid-r column density profile along z, plus max transverse (axis-0)
    # velocity. Keep `sim` alive by returning it.
    v = nodes.getFieldVector2d("velocity")
    imid = nr // 2
    prof = [density[grid.index(imid, j, 0)] for j in range(nz)]
    vtmax = max(abs(v[grid.index(i, j, 0)].x) for i in range(nr) for j in range(nz))
    return prof, vtmax, sim


def test_pure_z_reduces_to_cartesian():
    nr, nz, dr, dz, cycles = 12, 40, 1.0, 1.0, 150
    rz_prof, rz_vt, rz_sim = z_sod_profile(Geometry.CylindricalRZ, nr, nz, dr, dz, cycles)
    ca_prof, ca_vt, ca_sim = z_sod_profile(Geometry.Cartesian, nr, nz, dr, dz, cycles)

    # (i) the z-dynamics reduce to Cartesian: mid-column density profiles agree.
    maxdiff = max(abs(a - b) for a, b in zip(rz_prof, ca_prof))
    assert maxdiff < 1e-3, f"RZ z-Sod mid column differs from Cartesian by {maxdiff:.3e}"

    # (ii) no NEW radial instability: RZ radial noise stays within ~2x the base
    # scheme's transverse noise in Cartesian. Both sit at round-off once the
    # integrators refresh the ghost rim every sub-stage, so this needs an
    # absolute floor too -- a bare ratio is meaningless against a zero baseline.
    tol = max(2.0 * ca_vt, 1e-12)
    assert rz_vt < tol, (f"RZ radial noise {rz_vt:.3e} exceeds tolerance {tol:.3e} "
                         f"(2x Cartesian transverse noise {ca_vt:.3e})")
    print(f"[B ok] pure-z reduces to Cartesian: |RZ-Cart|_density={maxdiff:.2e}; "
          f"radial noise {rz_vt:.2e} vs Cartesian transverse {ca_vt:.2e} (tol {tol:.2e})")


def test_axis_guard_and_autoinstall():
    # Auto-install: a valid RZ run (no user 'left') just works and the axis
    # boundary is present (density stays bounded, no NaN at r=0).
    nr, nz = 16, 8
    sim = build(Geometry.CylindricalRZ, nr, nz, 0.05, 0.1)
    nodes = sim["nodes"]
    d = nodes.getFieldDouble("density"); e = nodes.getFieldDouble("specificInternalEnergy")
    for i in range(nr * nz):
        d.setValue(i, 1.0); e.setValue(i, 2.5)
    integ = integrator_for(sim["hydro"])
    Controller(integrator=integ, periodicWork=[], statStep=100000).Step(5)
    assert all(d[i] == d[i] for i in range(nr * nz)), "NaN near axis with auto-installed BC"
    print("[C ok] axis boundary auto-installed (no user 'left'): stable at r=0")

    # Guard: a user boundary that also claims 'left' must be rejected.
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
        print(f"[C ok] axis guard rejected user 'left' face: {str(ex)[:60]}...")
    assert raised, "expected RuntimeError when a user boundary claims the r=0 axis"


if __name__ == "__main__":
    test_uniform_rest()
    test_pure_z_reduces_to_cartesian()
    test_axis_guard_and_autoinstall()
    print("\nAll RZ invariants passed.")
