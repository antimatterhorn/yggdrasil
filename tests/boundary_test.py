import sys
from yggdrasil import *
from Mesh import Grid2d
from Physics import GridHydroKT2d
from EOS import IdealGasEOS
from Boundaries import (ReflectingGridBoundary2d,
                         OutflowGridBoundary2d,
                         PeriodicGridBoundary2d)

GAMMA = 1.4

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_hydro(nx, ny, dx=0.01, dy=0.01):
    # Return constants and eos explicitly so the caller keeps them alive.
    # The C++ hydro stores raw pointers to both; if they're GC'd early the
    # next physics call will segfault.
    grid      = Grid2d(nx, ny, dx, dy)
    nodeList  = NodeList(nx * ny)
    constants = MKS()
    eos       = IdealGasEOS(GAMMA, constants)
    hydro     = GridHydroKT2d(nodeList, constants, eos, grid)
    return grid, nodeList, hydro, constants, eos

def fill_uniform(nodeList, grid, nx, ny, rho=1.4, vx=1.0, vy=1.0, p=1.0):
    u  = p / ((GAMMA - 1.0) * rho)
    dn = nodeList.getFieldDouble("density")
    en = nodeList.getFieldDouble("specificInternalEnergy")
    vn = nodeList.getFieldVector2d("velocity")
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            dn.setValue(idx, rho)
            en.setValue(idx, u)
            vn.setValue(idx, Vector2d(vx, vy))

def fill_gradient(nodeList, grid, nx, ny, rho_fn, vx=0.0, vy=0.0, p=1.0):
    dn = nodeList.getFieldDouble("density")
    en = nodeList.getFieldDouble("specificInternalEnergy")
    vn = nodeList.getFieldVector2d("velocity")
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            rho = rho_fn(i, j)
            dn.setValue(idx, rho)
            en.setValue(idx, p / ((GAMMA - 1.0) * rho))
            vn.setValue(idx, Vector2d(vx, vy))

def run_one_step(hydro, dtmin=1e-6):
    ctrl = Controller(
        integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False),
        statStep   = 10000)
    ctrl.Step(1)

# ---------------------------------------------------------------------------
# Test runner
# ---------------------------------------------------------------------------

_results = []

def check(label, condition):
    status = "PASS" if condition else "FAIL"
    print(f"    [{status}] {label}")
    _results.append(condition)
    return condition

# ---------------------------------------------------------------------------
# Test 1: Reflecting — all faces (baseline)
# ---------------------------------------------------------------------------
def test_reflecting_all_faces():
    print("Test 1: Reflecting BC — all faces (default)")
    nx, ny = 7, 7
    grid, nodeList, hydro, *_ = make_hydro(nx, ny)
    fill_uniform(nodeList, grid, nx, ny, rho=1.4, vx=1.0, vy=1.0)
    bc = ReflectingGridBoundary2d(grid=grid)
    hydro.addBoundary(bc)
    run_one_step(hydro)
    v   = nodeList.getFieldVector2d("velocity")
    mid = 3
    check("left  ghost: vx flipped  (< 0)", v[grid.index(0,    mid, 0)].x < 0)
    check("left  ghost: vy intact   (> 0)", v[grid.index(0,    mid, 0)].y > 0)
    check("right ghost: vx flipped  (< 0)", v[grid.index(nx-1, mid, 0)].x < 0)
    check("right ghost: vy intact   (> 0)", v[grid.index(nx-1, mid, 0)].y > 0)
    check("bot   ghost: vy flipped  (< 0)", v[grid.index(mid,  0,   0)].y < 0)
    check("bot   ghost: vx intact   (> 0)", v[grid.index(mid,  0,   0)].x > 0)
    check("top   ghost: vy flipped  (< 0)", v[grid.index(mid, ny-1, 0)].y < 0)
    check("top   ghost: vx intact   (> 0)", v[grid.index(mid, ny-1, 0)].x > 0)

# ---------------------------------------------------------------------------
# Test 2: Reflecting — bottom face only via setFaces
# ---------------------------------------------------------------------------
def test_reflecting_bottom_only():
    print("Test 2: Reflecting BC — setFaces(['bottom'])")
    nx, ny = 7, 7
    grid, nodeList, hydro, *_ = make_hydro(nx, ny)
    fill_uniform(nodeList, grid, nx, ny, rho=1.4, vx=1.0, vy=1.0)
    bc = ReflectingGridBoundary2d(grid=grid)
    bc.setFaces(["bottom"])
    hydro.addBoundary(bc)
    run_one_step(hydro)
    v   = nodeList.getFieldVector2d("velocity")
    mid = 3
    check("bot   ghost: vy flipped  (< 0)", v[grid.index(mid,  0,   0)].y < 0)
    check("left  ghost: vx intact   (> 0)", v[grid.index(0,    mid, 0)].x > 0)
    check("left  ghost: vy intact   (> 0)", v[grid.index(0,    mid, 0)].y > 0)
    check("right ghost: vx intact   (> 0)", v[grid.index(nx-1, mid, 0)].x > 0)
    check("top   ghost: vy intact   (> 0)", v[grid.index(mid, ny-1, 0)].y > 0)

# ---------------------------------------------------------------------------
# Test 3: Reflecting — left face only via setFaces
# ---------------------------------------------------------------------------
def test_reflecting_left_only():
    print("Test 3: Reflecting BC — setFaces(['left'])")
    nx, ny = 7, 7
    grid, nodeList, hydro, *_ = make_hydro(nx, ny)
    fill_uniform(nodeList, grid, nx, ny, rho=1.4, vx=1.0, vy=1.0)
    bc = ReflectingGridBoundary2d(grid=grid)
    bc.setFaces(["left"])
    hydro.addBoundary(bc)
    run_one_step(hydro)
    v   = nodeList.getFieldVector2d("velocity")
    mid = 3
    check("left  ghost: vx flipped  (< 0)", v[grid.index(0,    mid, 0)].x < 0)
    check("right ghost: vx intact   (> 0)", v[grid.index(nx-1, mid, 0)].x > 0)
    check("bot   ghost: vy intact   (> 0)", v[grid.index(mid,  0,   0)].y > 0)
    check("top   ghost: vy intact   (> 0)", v[grid.index(mid, ny-1, 0)].y > 0)

# ---------------------------------------------------------------------------
# Test 4: Outflow — right face only via setFaces
#   x-only density gradient at uniform pressure → no interior pressure forces.
#   Outflow applies a Neumann (zero-gradient) copy: ghost ← adjacent interior.
#   Right ghost IC = 1.6, after BC ≈ 1.5 (pulled to match inner edge).
#   Left/top ghosts hold their IC values (no BC applied to them).
# ---------------------------------------------------------------------------
def test_outflow_right_only():
    print("Test 4: Outflow BC — setFaces(['right'])")
    nx, ny = 7, 7
    grid, nodeList, hydro, *_ = make_hydro(nx, ny)
    fill_gradient(nodeList, grid, nx, ny, rho_fn=lambda i, _: 1.0 + 0.1 * i)
    bc = OutflowGridBoundary2d(grid=grid)
    bc.setFaces(["right"])
    hydro.addBoundary(bc)
    run_one_step(hydro)
    d   = nodeList.getFieldDouble("density")
    mid = 3
    rho_ghost  = d[grid.index(nx-1, mid, 0)]
    rho_inner1 = d[grid.index(nx-2, mid, 0)]
    check("right ghost pulled from IC 1.6 toward inner ≈ 1.5", rho_ghost < 1.55)
    check("right ghost ≈ adjacent interior (Neumann copy)",
          abs(rho_ghost - rho_inner1) < 0.15)
    check("left  ghost: unchanged from IC (≈ 1.0)",
          abs(d[grid.index(0, mid, 0)] - 1.0) < 0.2)
    check("top   ghost: unchanged from IC",
          abs(d[grid.index(mid, ny-1, 0)] - (1.0 + 0.1 * mid)) < 0.2)

# ---------------------------------------------------------------------------
# Test 5: Periodic — setFaces(['left']) auto-pairs with right; top/bottom skip
# ---------------------------------------------------------------------------
def test_periodic_left_autopairs_right():
    print("Test 5: Periodic BC — setFaces(['left']) auto-pairs right, skips top/bottom")
    nx, ny = 7, 7
    grid, nodeList, hydro, *_ = make_hydro(nx, ny)
    rho_L, rho_R = 1.0, 3.0
    fill_gradient(nodeList, grid, nx, ny,
                  rho_fn=lambda i, j: rho_L if i < nx // 2 else rho_R)
    bc = PeriodicGridBoundary2d(grid=grid)
    bc.setFaces(["left"])
    hydro.addBoundary(bc)
    run_one_step(hydro)
    d   = nodeList.getFieldDouble("density")
    mid = 3   # i=3 is in the right half (rho_R), so j=mid ghost should ≈ rho_R
    check("left  ghost copies from right interior (> 2)",
          d[grid.index(0,    mid, 0)] > 2.0)
    check("right ghost copies from left interior  (< 2)",
          d[grid.index(nx-1, mid, 0)] < 2.0)
    check("top   ghost: not periodically wrapped (≈ rho_R)",
          d[grid.index(mid, ny-1, 0)] > 2.0)
    check("bot   ghost: not periodically wrapped (≈ rho_R)",
          d[grid.index(mid, 0,   0)] > 2.0)

# ---------------------------------------------------------------------------
# Test 6: Multiple BCs on different faces
# ---------------------------------------------------------------------------
def test_multi_bc():
    print("Test 6: Multi-BC — reflecting bottom + outflow left/right/top")
    nx, ny = 7, 7
    grid, nodeList, hydro, *_ = make_hydro(nx, ny)
    fill_uniform(nodeList, grid, nx, ny, rho=1.4, vx=1.0, vy=1.0)
    wall = ReflectingGridBoundary2d(grid=grid)
    wall.setFaces(["bottom"])
    hydro.addBoundary(wall)
    flow = OutflowGridBoundary2d(grid=grid)
    flow.setFaces(["left", "right", "top"])
    hydro.addBoundary(flow)
    run_one_step(hydro)
    v   = nodeList.getFieldVector2d("velocity")
    mid = 3
    check("bot   ghost: vy flipped  (< 0)", v[grid.index(mid,  0,   0)].y < 0)
    check("top   ghost: vy outflow  (> 0)", v[grid.index(mid, ny-1, 0)].y > 0)
    check("left  ghost: vx outflow  (> 0)", v[grid.index(0,    mid, 0)].x > 0)
    check("right ghost: vx outflow  (> 0)", v[grid.index(nx-1, mid, 0)].x > 0)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    test_reflecting_all_faces()
    test_reflecting_bottom_only()
    test_reflecting_left_only()
    test_outflow_right_only()
    test_periodic_left_autopairs_right()
    test_multi_bc()

    n_pass  = sum(_results)
    n_total = len(_results)
    print()
    if all(_results):
        print(f"ALL {n_total} CHECKS PASSED")
    else:
        print(f"{n_total - n_pass} / {n_total} CHECKS FAILED")
    sys.exit(0 if all(_results) else 1)
