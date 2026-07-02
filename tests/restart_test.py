import sys
import os
import tempfile
from yggdrasil import *
from Mesh import Grid2d
from Physics import GridHydroKT2d
from EOS import IdealGasEOS
from Boundaries import PeriodicGridBoundary2d
from IO import RestartWriter2d, RestartReader2d

GAMMA = 1.4

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_kh_problem(nx, ny, dx=0.02, dy=0.02, dtmin=1e-6):
    # Return every long-lived object explicitly -- the C++ side holds raw
    # pointers to constants/eos, so they must outlive the hydro package.
    grid      = Grid2d(nx, ny, dx, dy)
    nodeList  = NodeList(nx * ny)
    constants = MKS()
    eos       = IdealGasEOS(GAMMA, constants)
    hydro     = GridHydroKT2d(nodeList, constants, eos, grid)
    box       = PeriodicGridBoundary2d(grid=grid)
    hydro.addBoundary(box)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    density  = nodeList.getFieldDouble("density")
    energy   = nodeList.getFieldDouble("specificInternalEnergy")
    velocity = nodeList.getFieldVector2d("velocity")

    # Deterministic Kelvin-Helmholtz-style IC (no RNG, so two independent
    # constructions with the same nx/ny always agree exactly).
    p0 = 2.5
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            if j < ny // 4 or j >= 3 * ny // 4:
                rho, vx = 1.0, -0.5
            else:
                rho, vx = 2.0, 0.5
            vy = 0.1 if (i + j) % 2 == 0 else -0.1
            velocity.setValue(idx, Vector2d(vx, vy))
            density.setValue(idx, rho)
            energy.setValue(idx, p0 / ((GAMMA - 1.0) * rho))

    return grid, nodeList, hydro, integrator, constants, eos, box

def fields_match(nodeListA, nodeListB, names, tol=1e-10):
    for name in names:
        if name == "velocity":
            a = nodeListA.getFieldVector2d(name)
            b = nodeListB.getFieldVector2d(name)
            for i in range(a.size()):
                if abs(a[i].x - b[i].x) > tol or abs(a[i].y - b[i].y) > tol:
                    return False
        else:
            a = nodeListA.getFieldDouble(name)
            b = nodeListB.getFieldDouble(name)
            for i in range(a.size()):
                if abs(a[i] - b[i]) > tol:
                    return False
    return True

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
# Test 1: checkpoint + restore reproduces a continuous run exactly
# ---------------------------------------------------------------------------
def test_restart_round_trip():
    print("Test 1: checkpoint + restore reproduces a continuous run")
    nx, ny, cycles, half = 16, 16, 20, 10
    restartFile = os.path.join(tempfile.gettempdir(), "yggdrasil_restart_test.ygr")
    if os.path.exists(restartFile):
        os.remove(restartFile)

    # Ground truth: run straight through, dropping (but not using) a
    # checkpoint at the halfway point.
    grid, nodeList, hydro, integrator, *_keepalive = make_kh_problem(nx, ny)
    ctrl = Controller(integrator=integrator, statStep=10000)
    ctrl.Step(half)
    RestartWriter2d(nodeList, integrator).write(restartFile)
    ctrl.Step(cycles - half)
    truthCycle = integrator.Cycle()
    truthTime  = integrator.Time()

    # Restart run: fresh construction, restore from the halfway checkpoint,
    # then finish the remaining steps.
    grid2, nodeList2, hydro2, integrator2, *_keepalive2 = make_kh_problem(nx, ny)
    reader = RestartReader2d(nodeList2, integrator2)
    reader.read(restartFile)
    check("cycle restored to halfway point", integrator2.Cycle() == half)
    ctrl2 = Controller(integrator=integrator2, statStep=10000)
    ctrl2.Step(cycles - half)

    check("cycle matches continuous run", integrator2.Cycle() == truthCycle)
    check("time matches continuous run", abs(integrator2.Time() - truthTime) < 1e-12)
    check("density field matches continuous run",
          fields_match(nodeList, nodeList2, ["density"]))
    check("specificInternalEnergy field matches continuous run",
          fields_match(nodeList, nodeList2, ["specificInternalEnergy"]))
    check("velocity field matches continuous run",
          fields_match(nodeList, nodeList2, ["velocity"]))

    os.remove(restartFile)

# ---------------------------------------------------------------------------
# Test 2: mismatched NodeList size is rejected, not silently misapplied
# ---------------------------------------------------------------------------
def test_restart_rejects_size_mismatch():
    print("Test 2: restart rejects a NodeList of the wrong size")
    nx, ny = 8, 8
    restartFile = os.path.join(tempfile.gettempdir(), "yggdrasil_restart_test_mismatch.ygr")
    if os.path.exists(restartFile):
        os.remove(restartFile)

    grid, nodeList, hydro, integrator, *_keepalive = make_kh_problem(nx, ny)
    ctrl = Controller(integrator=integrator, statStep=10000)
    ctrl.Step(5)
    RestartWriter2d(nodeList, integrator).write(restartFile)

    grid2, nodeList2, hydro2, integrator2, *_keepalive2 = make_kh_problem(nx, ny + 2)
    raised = False
    try:
        RestartReader2d(nodeList2, integrator2).read(restartFile)
    except RuntimeError:
        raised = True
    check("mismatched NodeList size raises RuntimeError", raised)

    os.remove(restartFile)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    test_restart_round_trip()
    test_restart_rejects_size_mismatch()

    n_pass  = sum(_results)
    n_total = len(_results)
    print()
    if all(_results):
        print(f"ALL {n_total} CHECKS PASSED")
    else:
        print(f"{n_total - n_pass} / {n_total} CHECKS FAILED")
    sys.exit(0 if all(_results) else 1)
