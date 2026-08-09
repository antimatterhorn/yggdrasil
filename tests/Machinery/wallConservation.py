# Mass conservation at a domain boundary.
#
# Interior faces cancel exactly -- each is computed from identical arguments as
# both neighbours' flux_R/flux_L -- so any mass drift is the boundary's. Two
# independent mechanisms used to break that, and both are guarded here:
#
#   1. Ghost cells never receive a derivative, so an integrator's intermediate
#      state carries a rim one sub-stage stale. The integrators now call
#      Physics::ApplyStageBoundaries on every interim state.
#   2. At a domain-edge face only one side has a second-layer neighbour, so a
#      reconstructing solver took a one-sided slope there. Even with a perfectly
#      mirrored cell-centred pair, that made the reconstructed *face* states
#      asymmetric. GridHydroHLLC and GridHydroKT now drop both sides to
#      piecewise-constant at such a face; GridHydroHLLE never reconstructed and
#      so was only ever affected by (1).
#
# The periodic case is the same reconstruction bug seen from the other side: the
# two seam faces are one physical face, and a one-sided slope gave them
# different fluxes.

import sys
from yggdrasil import *
from Mesh import Grid2d
from EOS import IdealGasEOS
from Physics import GridHydroHLLC2d, GridHydroHLLE2d, GridHydroKT2d
from Boundaries import ReflectingGridBoundary2d, PeriodicGridBoundary2d

SOLVERS = [("HLLC", GridHydroHLLC2d),
           ("HLLE", GridHydroHLLE2d),
           ("KT",   GridHydroKT2d)]

GAMMA = 1.4
ROUNDOFF = 1e-13

_results = []

def check(label, condition):
    print(f"    [{'PASS' if condition else 'FAIL'}] {label}")
    _results.append(condition)
    return condition

def run(solver, integrator_cls, n, cycles, dx=0.05,
        boundary=ReflectingGridBoundary2d, offset=0.0):
    """Over-pressured blob relaxing into ambient gas inside a bounded domain.
    `offset` shifts the blob off-centre, in cells. Returns (relative mass drift
    over the real cells, max |vy| next to a wall)."""
    constants = MKS()
    eos       = IdealGasEOS(GAMMA, constants)
    grid      = Grid2d(n, n, dx, dx)
    nodeList  = NodeList(grid.size())
    hydro     = solver(nodeList, constants, eos, grid)
    hydro.addBoundary(boundary(grid=grid))
    integ = integrator_cls([hydro], dtmin=1e-8)

    cx = (0.5 * n + offset) * dx
    cy = (0.5 * n - offset) * dx
    rho = nodeList.getFieldDouble("density")
    u   = nodeList.getFieldDouble("specificInternalEnergy")
    v   = nodeList.getFieldVector2d("velocity")
    for i in range(grid.size()):
        pos  = grid.position(i)
        blob = (pos.x - cx) ** 2 + (pos.y - cy) ** 2 < (2.5 * dx) ** 2
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
# Test 1: every solver conserves mass while the wave is hitting the wall
# ---------------------------------------------------------------------------
def test_wall_conservation():
    print("Test 1: mass conserved at a reflecting wall being struck")
    for label, solver in SOLVERS:
        drift, wall = run(solver, RungeKutta2Integrator2d, 16, 40)
        check(f"{label} RK2: wave actually reaches the wall (|vy|={wall:.2e} > 1e-3)",
              wall > 1e-3)
        check(f"{label} RK2: mass drift {drift:+.2e} at round-off (< {ROUNDOFF:.0e})",
              abs(drift) < ROUNDOFF)

# ---------------------------------------------------------------------------
# Test 2: the stage-refresh half of the fix, which only RK sub-stages exercise
# ---------------------------------------------------------------------------
def test_wall_conservation_rk4():
    print("Test 2: same, under RK4 (three intermediate states per step)")
    for label, solver in SOLVERS:
        drift, wall = run(solver, RungeKutta4Integrator2d, 16, 40)
        check(f"{label} RK4: wave actually reaches the wall (|vy|={wall:.2e} > 1e-3)",
              wall > 1e-3)
        check(f"{label} RK4: mass drift {drift:+.2e} at round-off (< {ROUNDOFF:.0e})",
              abs(drift) < ROUNDOFF)

# ---------------------------------------------------------------------------
# Test 3: a periodic seam is one physical face and must not leak either
#
# The blob is offset so the solution is NOT symmetric about the domain centre.
# Centred, a periodic ghost and a reflected ghost carry identical values (the
# opposite-side cell *is* the mirror image), so the periodic case would silently
# reduce to Test 1 and prove nothing about the seam.
# ---------------------------------------------------------------------------
def test_periodic_seam_conservation():
    print("Test 3: mass conserved across a periodic seam (off-centre blob)")
    for label, solver in SOLVERS:
        drift, _ = run(solver, RungeKutta2Integrator2d, 16, 40,
                       boundary=PeriodicGridBoundary2d, offset=2.5)
        check(f"{label}: mass drift {drift:+.2e} at round-off (< {ROUNDOFF:.0e})",
              abs(drift) < ROUNDOFF)

# ---------------------------------------------------------------------------
# Test 4: control — with the wall untouched, every solver is exactly
# conservative, confirming a drift above could only come from the boundary
# ---------------------------------------------------------------------------
def test_interior_is_conservative():
    print("Test 4: control — interior-only flow conserves mass in all solvers")
    for label, solver in SOLVERS:
        drift, wall = run(solver, RungeKutta2Integrator2d, 48, 40)
        check(f"{label}: wall untouched (|vy|={wall:.2e})", wall < 1e-12)
        check(f"{label}: mass drift {drift:+.2e} at round-off (< {ROUNDOFF:.0e})",
              abs(drift) < ROUNDOFF)

if __name__ == "__main__":
    test_wall_conservation()
    test_wall_conservation_rk4()
    test_periodic_seam_conservation()
    test_interior_is_conservative()

    n_pass, n_total = sum(_results), len(_results)
    print()
    if all(_results):
        print(f"ALL {n_total} CHECKS PASSED")
    else:
        print(f"{n_total - n_pass} / {n_total} CHECKS FAILED")
    sys.exit(0 if all(_results) else 1)
