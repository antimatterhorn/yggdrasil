from yggdrasil import *
from Mesh import Geometry

# Cylindrical (r,z) Sedov point blast at the (r=0, z=0) corner of an
# axisymmetric grid, exercising the CylindricalRZ metric, the p/r source term,
# and the auto-installed r=0 axis boundary. Checks: the solution stays finite,
# the shock expands outward in radius, and per-radian mass (sum rho*cellVolume)
# is conserved -- reflecting walls + axis => no net flux out of the interior.

GAMMA = 1.4


def _mass(grid, density, ids):
    return sum(density[idx] * grid.cellVolume(idx) for idx in ids)


def run():
    nr, nz, dr = 60, 60, 0.02
    grid = Grid2d(nr, nz, dr, dr, Geometry.CylindricalRZ)
    nodes = NodeList(nr * nz)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroHLLE2d(nodes, constants, eos, grid)
    box = ReflectingGridBoundary2d(grid=grid)
    box.setFaces(["right", "top", "bottom"])          # r=0 axis auto-installed
    hydro.addBoundary(box)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=1e-6, verbose=False)

    amb = 1.0
    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    for i in range(nr * nz):
        density[i] = amb; energy[i] = 1e-3
    for j in range(2):
        for i in range(2):
            energy[grid.index(i, j, 0)] = 400.0

    ids = [grid.index(i, j, 0) for j in range(nz) for i in range(nr)
           if not grid.onBoundary(grid.index(i, j, 0))]

    def shock_radius():
        return max((grid.cellRadius(grid.index(i, j, 0))
                    for j in range(nz) for i in range(nr)
                    if density[grid.index(i, j, 0)] > 1.5 * amb), default=0.0)

    m0 = _mass(grid, density, ids)
    Controller(integrator=integrator, periodicWork=[], statStep=100000).Step(400)

    finite = all(density[i] == density[i] and abs(density[i]) < 1e6 for i in range(nr * nz))
    r1 = shock_radius()
    dm = abs(_mass(grid, density, ids) - m0) / m0

    return {"mode": "analytic", "checks": [
        ("solution finite and bounded", bool(finite), "no NaN/inf, rho < 1e6"),
        ("shock expanded radially", bool(r1 > 5 * dr), f"R={r1:.3f} > {5*dr:.3f}"),
        ("per-radian mass conserved", bool(dm < 1e-3), f"dM/M={dm:.2e} < 1e-3"),
    ]}


if __name__ == "__main__":
    print(run())
