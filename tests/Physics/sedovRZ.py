from yggdrasil import *
from Physics import GridHydroHLLE2d
from Mesh import Grid2d, Geometry
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d

# Sedov-like point blast at the (r=0, z=0) corner of an axisymmetric (r,z) grid.
# A concentrated pot of internal energy drives a shock that must expand into the
# ambient gas and remain stable at the r=0 axis. This exercises the full RZ path
# (radius-weighted volumes/areas, the p/r source, and the auto-installed axis
# boundary) on a genuinely radial problem. Checks:
#   - the solution stays finite (no NaN/inf, bounded density);
#   - the shock front moves outward in radius;
#   - per-radian mass (sum rho*cellVolume) is conserved to a small tolerance
#     (reflecting walls + axis => no net flux out of the interior).

GAMMA = 1.4


def interior_ids(grid, nr, nz):
    return [grid.index(i, j, 0) for j in range(nz) for i in range(nr)
            if not grid.onBoundary(grid.index(i, j, 0))]


def total_mass(density, grid, ids):
    return sum(density[idx] * grid.cellVolume(idx) for idx in ids)


def shock_radius(density, grid, nr, nz, rho_ambient):
    # outermost cell radius whose density exceeds 1.5x ambient (the shell)
    rmax = 0.0
    for j in range(nz):
        for i in range(nr):
            idx = grid.index(i, j, 0)
            if density[idx] > 1.5 * rho_ambient:
                rmax = max(rmax, grid.cellRadius(idx))
    return rmax


if __name__ == "__main__":
    nr = nz = 60
    dr = dz = 0.02
    grid = Grid2d(nr, nz, dr, dz, Geometry.CylindricalRZ)
    nodes = NodeList(nr * nz)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroHLLE2d(nodes, constants, eos, grid)
    box = ReflectingGridBoundary2d(grid=grid)
    box.setFaces(["right", "top", "bottom"])   # 'left' (r=0 axis) auto-installed
    hydro.addBoundary(box)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=1e-6, verbose=False)

    rho_ambient = 1.0
    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    for i in range(nr * nz):
        density.setValue(i, rho_ambient)
        energy.setValue(i, 1e-3)               # cold ambient

    # Deposit energy in a 2x2 blob at the axis/base corner.
    blob = 2
    E_blob = 400.0
    for j in range(blob):
        for i in range(blob):
            energy.setValue(grid.index(i, j, 0), E_blob)

    ids = interior_ids(grid, nr, nz)
    m0 = total_mass(density, grid, ids)
    r0 = shock_radius(density, grid, nr, nz, rho_ambient)

    controller = Controller(integrator=integrator, periodicWork=[], statStep=100000)
    controller.Step(400)

    # (i) finite and bounded
    finite = all(density[i] == density[i] and abs(density[i]) < 1e6
                 for i in range(nr * nz))
    assert finite, "blast produced NaN/inf or runaway density"
    print("[ok] solution finite and bounded after 400 steps")

    # (ii) shock expanded outward
    r1 = shock_radius(density, grid, nr, nz, rho_ambient)
    assert r1 > r0 + 5 * dr, f"shock did not expand: r0={r0:.3f} -> r1={r1:.3f}"
    print(f"[ok] shock front expanded radially: r={r0:.3f} -> {r1:.3f}")

    # (iii) per-radian mass conserved
    m1 = total_mass(density, grid, ids)
    rel = abs(m1 - m0) / m0
    assert rel < 1e-3, f"mass not conserved: rel change {rel:.3e}"
    print(f"[ok] per-radian mass conserved: rel change {rel:.2e}")

    print("\nCylindrical Sedov blast: all checks passed.")
