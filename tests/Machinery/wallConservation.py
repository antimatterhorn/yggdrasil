# Mass conservation at a reflecting wall.
#
# Ghost cells never receive a derivative, so an integrator's intermediate state
# carries a rim one sub-stage stale unless the boundaries are reapplied. At a
# reflecting wall that breaks the mirror symmetry the wall's Riemann problem
# depends on, and the domain leaks mass through it -- interior faces cancel
# exactly, so any drift is the wall's. The integrators fix this by calling
# Physics::ApplyStageBoundaries on every interim state.
#
# HLLE pairs the ghost against its interior partner with no slope
# reconstruction, so its wall flux is exactly zero once the rim is fresh and its
# drift must sit at round-off. HLLC and KT reconstruct a one-sided slope at a
# domain-edge face, which breaks the same symmetry independently of rim
# staleness; they are checked only on the control below, where the wave never
# reaches the wall.

import sys
from yggdrasil import *
from Mesh import Grid2d
from EOS import IdealGasEOS
from Physics import GridHydroHLLC2d, GridHydroHLLE2d, GridHydroKT2d
from Boundaries import ReflectingGridBoundary2d

GAMMA = 1.4
ROUNDOFF = 1e-13

_results = []

def check(label, condition):
    print(f"    [{'PASS' if condition else 'FAIL'}] {label}")
    _results.append(condition)
    return condition

def run(solver, integrator_cls, n, cycles, dx=0.05):
    """Over-pressured blob relaxing into ambient gas inside reflecting walls.
    Returns (relative mass drift over the real cells, max |vy| next to a wall)."""
    constants = MKS()
    eos       = IdealGasEOS(GAMMA, constants)
    grid      = Grid2d(n, n, dx, dx)
    nodeList  = NodeList(grid.size())
    hydro     = solver(nodeList, constants, eos, grid)
    hydro.addBoundary(ReflectingGridBoundary2d(grid=grid))
    integ = integrator_cls([hydro], dtmin=1e-8)

    centre = 0.5 * n * dx
    rho = nodeList.getFieldDouble("density")
    u   = nodeList.getFieldDouble("specificInternalEnergy")
    v   = nodeList.getFieldVector2d("velocity")
    for i in range(grid.size()):
        pos  = grid.position(i)
        blob = (pos.x - centre) ** 2 + (pos.y - centre) ** 2 < (2.5 * dx) ** 2
        p, r = (5.0, 2.0) if blob else (1.0, 1.0)
        rho.setValue(i, r)
        u.setValue(i, p / ((GAMMA - 1.0) * r))
        v.setValue(i, Vector2d(0.0, 0.0))

    real = [grid.index(i, j, 0) for j in range(1, n - 1) for i in range(1, n - 1)]
    mass = lambda: sum(rho[i] * grid.cellVolume(i) for i in real)

    integ.Initialize()
    m0 = mass()
    for _ in range(cycles):
        integ.Step()
    wall = max(abs(v[grid.index(i, 1, 0)].y) for i in range(1, n - 1))
    return (mass() - m0) / m0, wall

# ---------------------------------------------------------------------------
# Test 1: HLLE conserves mass exactly while the wave is hitting the wall
# ---------------------------------------------------------------------------
def test_hlle_wall_conservation():
    print("Test 1: HLLE — mass conserved at a reflecting wall being struck")
    for name, integrator_cls in [("RK2", RungeKutta2Integrator2d),
                                 ("RK4", RungeKutta4Integrator2d)]:
        drift, wall = run(GridHydroHLLE2d, integrator_cls, 16, 40)
        check(f"{name}: wave actually reaches the wall (|vy|={wall:.2e} > 1e-3)",
              wall > 1e-3)
        check(f"{name}: mass drift {drift:+.2e} at round-off (< {ROUNDOFF:.0e})",
              abs(drift) < ROUNDOFF)

# ---------------------------------------------------------------------------
# Test 2: control — with the wall untouched, every solver is exactly
# conservative, confirming a drift in Test 1 could only come from the wall
# ---------------------------------------------------------------------------
def test_interior_is_conservative():
    print("Test 2: control — interior-only flow conserves mass in all solvers")
    for label, solver in [("HLLC", GridHydroHLLC2d),
                          ("HLLE", GridHydroHLLE2d),
                          ("KT",   GridHydroKT2d)]:
        drift, wall = run(solver, RungeKutta2Integrator2d, 48, 40)
        check(f"{label}: wall untouched (|vy|={wall:.2e})", wall < 1e-12)
        check(f"{label}: mass drift {drift:+.2e} at round-off (< {ROUNDOFF:.0e})",
              abs(drift) < ROUNDOFF)

if __name__ == "__main__":
    test_hlle_wall_conservation()
    test_interior_is_conservative()

    n_pass, n_total = sum(_results), len(_results)
    print()
    if all(_results):
        print(f"ALL {n_total} CHECKS PASSED")
    else:
        print(f"{n_total - n_pass} / {n_total} CHECKS FAILED")
    sys.exit(0 if all(_results) else 1)
