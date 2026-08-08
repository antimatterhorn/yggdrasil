import sys
from yggdrasil import *
from Mesh import Grid2d, Grid3d
from Physics import GridHydroKT2d, GridHydroKT3d
from EOS import IdealGasEOS
from Boundaries import (ReflectingGridBoundary2d,
                         OutflowGridBoundary2d,
                         PeriodicGridBoundary2d,
                         PeriodicGridBoundary3d,
                         OutflowGridBoundary3d,
                         ReflectingGridBoundary3d,
                         DirichletGridBoundary3d)

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

def make_hydro_3d(nx, ny, nz, dx=0.01, dy=0.01, dz=0.01):
    grid      = Grid3d(nx, ny, nz, dx, dy, dz)
    nodeList  = NodeList(nx * ny * nz)
    constants = MKS()
    eos       = IdealGasEOS(GAMMA, constants)
    hydro     = GridHydroKT3d(nodeList, constants, eos, grid)
    return grid, nodeList, hydro, constants, eos

def fill_gradient_3d(nodeList, grid, nx, ny, nz, rho_fn, p=1.0):
    dn = nodeList.getFieldDouble("density")
    en = nodeList.getFieldDouble("specificInternalEnergy")
    vn = nodeList.getFieldVector3d("velocity")
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                idx = grid.index(i, j, k)
                rho = rho_fn(i, j, k)
                dn.setValue(idx, rho)
                en.setValue(idx, p / ((GAMMA - 1.0) * rho))
                vn.setValue(idx, Vector3d(0.0, 0.0, 0.0))

def run_one_step_3d(hydro, dtmin=1e-6):
    ctrl = Controller(
        integrator = RungeKutta4Integrator3d([hydro], dtmin=dtmin, verbose=False),
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
# Test 7: Periodic — 3D, non-cubic grid, all three axis pairs
#
# Regression test for two bugs that both required a non-cubic 3D grid to
# surface (a cubic grid masks them since nx==ny==nz numerically):
#   1. Grid<3>::findBoundaries computed the high-x face list using size()-1-b
#      (total cell count) instead of nx-1-b, producing wildly out-of-range
#      indices that PeriodicGridBoundary<3> fed straight into an unchecked
#      Field::setValue -- a heap buffer overflow write.
#   2. PeriodicGridBoundary<3>'s constructor never read grid->size_z() and
#      reused Nx/Ny for loop bounds that should have involved Nz, silently
#      mismatching the top/bottom and front/back ghost-copy pairings.
# ---------------------------------------------------------------------------
def test_periodic_3d_noncubic():
    print("Test 7: Periodic BC — 3D, non-cubic grid (all axes)")
    nx, ny, nz = 8, 6, 4
    rho_L, rho_R = 1.0, 3.0
    mid_i, mid_j, mid_k = nx // 2, ny // 2, nz // 2

    # x-axis: left (i=0) / right (i=nx-1)
    grid, nodeList, hydro, *_ = make_hydro_3d(nx, ny, nz)
    fill_gradient_3d(nodeList, grid, nx, ny, nz,
                      rho_fn=lambda i, j, k: rho_L if i < mid_i else rho_R)
    hydro.addBoundary(PeriodicGridBoundary3d(grid=grid))
    run_one_step_3d(hydro)
    d = nodeList.getFieldDouble("density")
    check("x: i=0 ghost copies from i=nx-1 side interior (> 2)",
          d[grid.index(0, mid_j, mid_k)] > 2.0)
    check("x: i=nx-1 ghost copies from i=0 side interior (< 2)",
          d[grid.index(nx - 1, mid_j, mid_k)] < 2.0)

    # y-axis: j=0 / j=ny-1
    grid, nodeList, hydro, *_ = make_hydro_3d(nx, ny, nz)
    fill_gradient_3d(nodeList, grid, nx, ny, nz,
                      rho_fn=lambda i, j, k: rho_L if j < mid_j else rho_R)
    hydro.addBoundary(PeriodicGridBoundary3d(grid=grid))
    run_one_step_3d(hydro)
    d = nodeList.getFieldDouble("density")
    check("y: j=0 ghost copies from j=ny-1 side interior (> 2)",
          d[grid.index(mid_i, 0, mid_k)] > 2.0)
    check("y: j=ny-1 ghost copies from j=0 side interior (< 2)",
          d[grid.index(mid_i, ny - 1, mid_k)] < 2.0)

    # z-axis: k=0 / k=nz-1
    grid, nodeList, hydro, *_ = make_hydro_3d(nx, ny, nz)
    fill_gradient_3d(nodeList, grid, nx, ny, nz,
                      rho_fn=lambda i, j, k: rho_L if k < mid_k else rho_R)
    hydro.addBoundary(PeriodicGridBoundary3d(grid=grid))
    run_one_step_3d(hydro)
    d = nodeList.getFieldDouble("density")
    check("z: k=0 ghost copies from k=nz-1 side interior (> 2)",
          d[grid.index(mid_i, mid_j, 0)] > 2.0)
    check("z: k=nz-1 ghost copies from k=0 side interior (< 2)",
          d[grid.index(mid_i, mid_j, nz - 1)] < 2.0)

# ---------------------------------------------------------------------------
# Test 8: Outflow — 3D, non-cubic grid, all three axis pairs
#
# Regression test for the same class of bug as Test 7, independently
# duplicated in OutflowGridBoundary<3>'s constructor: missing grid->size_z()
# and Nx/Ny reused for loop bounds that should have involved Nz.
# ---------------------------------------------------------------------------
def test_outflow_3d_noncubic():
    print("Test 8: Outflow BC — 3D, non-cubic grid (all axes)")
    nx, ny, nz = 8, 6, 4
    mid_i, mid_j, mid_k = nx // 2, ny // 2, nz // 2

    # Density varies only along x; uniform along y/z. A correctly axis-paired
    # outflow BC should pull the x-boundary ghosts toward the interior
    # x-gradient, while leaving y/z-boundary ghosts close to their (uniform)
    # IC -- if front/back or top/bottom accidentally read along the wrong
    # axis, these would pick up the x-gradient instead.
    grid, nodeList, hydro, *_ = make_hydro_3d(nx, ny, nz)
    fill_gradient_3d(nodeList, grid, nx, ny, nz, rho_fn=lambda i, j, k: 1.0 + 0.1 * i)
    hydro.addBoundary(OutflowGridBoundary3d(grid=grid))
    run_one_step_3d(hydro)
    d = nodeList.getFieldDouble("density")

    rho_ghost_right = d[grid.index(nx - 1, mid_j, mid_k)]
    rho_inner_right = d[grid.index(nx - 2, mid_j, mid_k)]
    check("x: right ghost pulled toward interior x-gradient",
          abs(rho_ghost_right - rho_inner_right) < 0.3)

    ic_uniform = 1.0 + 0.1 * mid_i
    check("y: j=0 ghost stays close to uniform-along-y IC (no cross-axis mixup)",
          abs(d[grid.index(mid_i, 0, mid_k)] - ic_uniform) < 0.3)
    check("z: k=0 ghost stays close to uniform-along-z IC (no cross-axis mixup)",
          abs(d[grid.index(mid_i, mid_j, 0)] - ic_uniform) < 0.3)

# ---------------------------------------------------------------------------
# Test 9: Reflecting — 3D, non-cubic grid, correct per-cell neighbor pairing
#
# Regression test for a third variant of the same underlying bug, in
# ReflectingGridBoundary<3>'s constructor: Nx/Ny/Nz were all correctly read
# via size_x()/size_y()/size_z(), but the interior-neighbor loops were nested
# in the wrong order (e.g. k outer / j inner instead of Grid<3>::
# findBoundaries' j outer / k inner for left-right). Position i in
# boundaryLists and position i in interiorLists must refer to the *same*
# (j,k)/(i,k)/(i,j) cell -- with mismatched nesting this only coincides when
# the two axis extents happen to be equal, so a non-cubic grid (nx=8,ny=6,
# nz=4, all different) is required to catch it.
#
# velocity is seeded as (1.0, 10*j, 10*k) so the flipped (normal) component
# is uninformative, but the two *preserved* components directly encode the
# (j,k)/(i,k)/(i,j) position of whichever interior cell actually got copied --
# a wrong-order pairing would show up as a value tied to some other cell's
# coordinates instead of the ghost's own.
# ---------------------------------------------------------------------------
def test_reflecting_3d_noncubic_pairing():
    print("Test 9: Reflecting BC — 3D, non-cubic grid, correct neighbor pairing")
    nx, ny, nz = 8, 6, 4
    grid, nodeList, hydro, *_ = make_hydro_3d(nx, ny, nz)
    velocity = nodeList.getFieldVector3d("velocity")
    density  = nodeList.getFieldDouble("density")
    energy   = nodeList.getFieldDouble("specificInternalEnergy")
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                idx = grid.index(i, j, k)
                velocity.setValue(idx, Vector3d(1.0, 10.0 * j, 10.0 * k))
                density.setValue(idx, 1.0)
                energy.setValue(idx, 2.5 / (GAMMA - 1.0))
    hydro.addBoundary(ReflectingGridBoundary3d(grid=grid))
    # dtmin deliberately tiny: keeps hydro-driven evolution of the seeded
    # velocity gradient negligible next to the coefficient (10) used above,
    # so any leftover drift can't be confused with a genuine mis-pairing.
    ctrl = Controller(
        integrator = RungeKutta4Integrator3d([hydro], dtmin=1e-9, verbose=False),
        statStep   = 10000)
    ctrl.Step(1)
    v = nodeList.getFieldVector3d("velocity")

    j0, k0 = 2, 1
    vL = v[grid.index(0, j0, k0)]
    check("left  ghost preserves its own (j,k) -- vy",  abs(vL.y - 10.0 * j0) < 1.0)
    check("left  ghost preserves its own (j,k) -- vz",  abs(vL.z - 10.0 * k0) < 1.0)

    i0, k1 = 3, 2
    vT = v[grid.index(i0, 0, k1)]
    check("top   ghost preserves its own (i,k) -- vx",  abs(vT.x - 1.0) < 1.0)
    check("top   ghost preserves its own (i,k) -- vz",  abs(vT.z - 10.0 * k1) < 1.0)

    i1, j1 = 4, 3
    vF = v[grid.index(i1, j1, 0)]
    check("front ghost preserves its own (i,j) -- vx",  abs(vF.x - 1.0) < 1.0)
    check("front ghost preserves its own (i,j) -- vy",  abs(vF.y - 10.0 * j1) < 1.0)

# ---------------------------------------------------------------------------
# Test 10: Dirichlet — 3D, non-cubic grid
#
# DirichletGridBoundary<3> never builds a separately-ordered interior-neighbor
# list (it just saves/restores IC values at the same `ids` list from
# lowMost()/highMost()), so it isn't exposed to the Nx/Ny/Nz-mislabeling
# bug the other three had. It WAS still transitively exposed to the root-cause
# bug in Grid<3>::findBoundaries() (the high-x face list returning out-of-range
# indices on any 3D grid, cubic or not) via its own use of it, so it's worth
# its own non-cubic regression check now that that's fixed.
# ---------------------------------------------------------------------------
def test_dirichlet_3d_noncubic():
    print("Test 10: Dirichlet BC — 3D, non-cubic grid")
    nx, ny, nz = 8, 6, 4
    grid, nodeList, hydro, *_ = make_hydro_3d(nx, ny, nz)
    fill_gradient_3d(nodeList, grid, nx, ny, nz, rho_fn=lambda i, j, k: 1.0 + 0.1 * i)
    bc = DirichletGridBoundary3d(grid=grid)
    bc.addDomain()
    check("Dirichlet ids all within grid bounds",
          max(bc.boundaryIds()) < grid.size())
    hydro.addBoundary(bc)
    run_one_step_3d(hydro)
    d = nodeList.getFieldDouble("density")
    mid_j, mid_k = ny // 2, nz // 2
    check("left  ghost pinned to IC (1.0)",
          abs(d[grid.index(0, mid_j, mid_k)] - 1.0) < 1e-6)
    check("right ghost pinned to IC (1.0 + 0.1*(nx-1))",
          abs(d[grid.index(nx - 1, mid_j, mid_k)] - (1.0 + 0.1 * (nx - 1))) < 1e-6)

# ---------------------------------------------------------------------------
# Test 11: Grid face lists are on the geometric side their axis argument names
#
# Grid used to expose these as leftMost()/rightMost()/topMost()/bottomMost(),
# where topMost() returned the *low*-j row -- inverted relative to how every
# caller and every setFaces() name uses "top". GridBoundary::setFaces indexes
# its names array [-x,+x,-y,+y,-z,+z], so the two inversions cancelled and the
# behaviour was right, but reading either half alone made it look like top and
# bottom were swapped. Pin the orientation down directly.
# ---------------------------------------------------------------------------
def test_grid_face_list_orientation():
    print("Test 11: Grid.lowMost/highMost are on the named side of each axis")
    nx, ny = 5, 8
    grid = Grid2d(nx, ny, 0.01, 0.01)
    check("lowMost(0)  is the i=0 column",
          sorted(grid.lowMost(0))  == sorted(grid.index(0, j, 0)      for j in range(ny)))
    check("highMost(0) is the i=nx-1 column",
          sorted(grid.highMost(0)) == sorted(grid.index(nx - 1, j, 0) for j in range(ny)))
    check("lowMost(1)  is the j=0 row",
          sorted(grid.lowMost(1))  == sorted(grid.index(i, 0, 0)      for i in range(nx)))
    check("highMost(1) is the j=ny-1 row",
          sorted(grid.highMost(1)) == sorted(grid.index(i, ny - 1, 0) for i in range(nx)))

# ---------------------------------------------------------------------------
# Test 12: every setFaces name activates the edge it is named after
#
# Only 'bottom' and 'left' were covered in isolation before (Tests 2 and 3);
# 'top' and 'right' were only ever exercised together with their partner face,
# which is exactly the case that hides a low/high mix-up.
# ---------------------------------------------------------------------------
def test_reflecting_each_face_orientation():
    print("Test 12: Reflecting BC — each setFaces name hits its own edge")
    nx, ny = 7, 7
    mid = 3

    # (face name, ghost cell that must flip, cell that must stay intact, component)
    cases = [
        ("left",   (0, mid),      (nx - 1, mid), "x"),
        ("right",  (nx - 1, mid), (0, mid),      "x"),
        ("bottom", (mid, 0),      (mid, ny - 1), "y"),
        ("top",    (mid, ny - 1), (mid, 0),      "y"),
    ]

    for face, flipped, intact, comp in cases:
        grid, nodeList, hydro, *_ = make_hydro(nx, ny)
        fill_uniform(nodeList, grid, nx, ny, rho=1.4, vx=1.0, vy=1.0)
        bc = ReflectingGridBoundary2d(grid=grid)
        bc.setFaces([face])
        hydro.addBoundary(bc)
        run_one_step(hydro)
        v = nodeList.getFieldVector2d("velocity")
        get = (lambda p: v[grid.index(p[0], p[1], 0)].x) if comp == "x" else \
              (lambda p: v[grid.index(p[0], p[1], 0)].y)
        check(f"setFaces(['{face}']): v{comp} flipped at {flipped}", get(flipped) < 0)
        check(f"setFaces(['{face}']): v{comp} intact  at {intact}",  get(intact)  > 0)

# ---------------------------------------------------------------------------
# Test 13: outflow top-only / bottom-only pick opposite rows
# ---------------------------------------------------------------------------
def test_outflow_top_and_bottom_only():
    print("Test 13: Outflow BC — setFaces(['top']) vs setFaces(['bottom'])")
    nx, ny = 7, 7
    mid = 3
    rho_fn = lambda i, j: 1.0 + 0.1 * j

    for face, ghost_j, inner_j, untouched_j in [("top",    ny - 1, ny - 2, 0),
                                                ("bottom", 0,      1,      ny - 1)]:
        grid, nodeList, hydro, *_ = make_hydro(nx, ny)
        fill_gradient(nodeList, grid, nx, ny, rho_fn)
        bc = OutflowGridBoundary2d(grid=grid)
        bc.setFaces([face])
        hydro.addBoundary(bc)
        run_one_step(hydro)
        d = nodeList.getFieldDouble("density")
        check(f"setFaces(['{face}']): j={ghost_j} ghost copies j={inner_j}",
              abs(d[grid.index(mid, ghost_j, 0)] - d[grid.index(mid, inner_j, 0)]) < 1e-12)
        check(f"setFaces(['{face}']): j={untouched_j} ghost left at its IC",
              abs(d[grid.index(mid, untouched_j, 0)] - rho_fn(mid, untouched_j)) < 1e-12)

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
    test_periodic_3d_noncubic()
    test_outflow_3d_noncubic()
    test_reflecting_3d_noncubic_pairing()
    test_dirichlet_3d_noncubic()
    test_grid_face_list_orientation()
    test_reflecting_each_face_orientation()
    test_outflow_top_and_bottom_only()

    n_pass  = sum(_results)
    n_total = len(_results)
    print()
    if all(_results):
        print(f"ALL {n_total} CHECKS PASSED")
    else:
        print(f"{n_total - n_pass} / {n_total} CHECKS FAILED")
    sys.exit(0 if all(_results) else 1)
