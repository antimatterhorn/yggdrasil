from yggdrasil import *
from Mesh import Grid2d
from EOS import IdealGasEOS
from Physics import GridHydroHLLC2d
from Boundaries import ReflectingGridBoundary2d
from AMRController import (AMRController, buildPatch, wireSiblingBoundaries,
                            wireCoarseFineBoundary, buildRestriction, buildReflux)

# Milestone 1: Patch/AMRController scaffolding. Checks that a patch built
# from an arbitrary (level, box) lands at the correct world-space position.

GAMMA = 1.4
REFINEMENT_RATIO = 2
NGHOST = 1

baseNx, baseNy = 8, 8
baseDx, baseDy = 1.0, 1.0

constants = MKS()
eos = IdealGasEOS(GAMMA, constants)

def physicsFactory(nodeList, grid):
    return GridHydroHLLC2d(nodeList, constants, eos, grid)

baseGrid = Grid2d(baseNx, baseNy, baseDx, baseDy)
baseNodeList = NodeList(baseGrid.size())

amr = AMRController(baseGrid, baseNodeList, physicsFactory,
                     RungeKutta2Integrator2d, dtmin=1e-7,
                     refinementRatio=REFINEMENT_RATIO, nghost=NGHOST)

# --- Level 0 is exactly the user's own grid/nodeList, untouched. ---
basePatch = amr.levels[0][0]
assert basePatch.grid is baseGrid
assert basePatch.nodeList is baseNodeList
assert basePatch.box == (0, 0, baseNx - 1, baseNy - 1)
print("[ok] level 0 patch wraps the user's base grid/nodeList unchanged.")

# --- Two same-level adjacent patches: verify contiguous world placement. ---
patchA = amr.buildPatch(level=1, box=(0, 0, 3, 7))
patchB = amr.buildPatch(level=1, box=(4, 0, 7, 7))

dx1, dy1 = amr.spacing(1)
assert abs(dx1 - baseDx / REFINEMENT_RATIO) < 1e-14
assert abs(dy1 - baseDy / REFINEMENT_RATIO) < 1e-14
assert patchA.grid.nx == (3 - 0 + 1) + 2 * NGHOST
assert patchA.grid.ny == (7 - 0 + 1) + 2 * NGHOST

# Rightmost interior cell of A and leftmost interior cell of B should be
# exactly one dx1 apart, center to center -- confirming no gap/overlap.
rightmostA = patchA.grid.index(NGHOST + 3, NGHOST + 0, 0)   # interior i=3 (box iHigh)
leftmostB  = patchB.grid.index(NGHOST + 0, NGHOST + 0, 0)   # interior i=0 (box iLow)
xA = patchA.grid.position(rightmostA).x
xB = patchB.grid.position(leftmostB).x
assert abs((xB - xA) - dx1) < 1e-12, f"patches not contiguous: xA={xA}, xB={xB}, dx1={dx1}"
print("[ok] two same-level patches built from adjacent boxes are contiguous in world space.")

# --- Parent/child nesting: fine cells must average exactly to the coarse
# cell center they refine, confirming the setOrigin sign convention. ---
ci, cj = 2, 3
coarseCenter = baseGrid.position(baseGrid.index(ci, cj, 0))

fineBox = (ci * REFINEMENT_RATIO, cj * REFINEMENT_RATIO,
           ci * REFINEMENT_RATIO + REFINEMENT_RATIO - 1,
           cj * REFINEMENT_RATIO + REFINEMENT_RATIO - 1)
finePatch = amr.buildPatch(level=1, box=fineBox)

fineCenters = []
for fj in range(REFINEMENT_RATIO):
    for fi in range(REFINEMENT_RATIO):
        idx = finePatch.grid.index(NGHOST + fi, NGHOST + fj, 0)
        fineCenters.append(finePatch.grid.position(idx))

avgX = sum(p.x for p in fineCenters) / len(fineCenters)
avgY = sum(p.y for p in fineCenters) / len(fineCenters)
assert abs(avgX - coarseCenter.x) < 1e-12, f"fine avg x={avgX} != coarse x={coarseCenter.x}"
assert abs(avgY - coarseCenter.y) < 1e-12, f"fine avg y={avgY} != coarse y={coarseCenter.y}"
print("[ok] a refinement-ratio^2 block of fine cells averages exactly to its parent coarse cell center.")

print("[ok] AMR milestone 1 (Patch/AMRController scaffolding) passes.")

# ============================================================================
# Milestone 2: AMR::PatchNeighborBoundary. Splits one domain into two patches
# joined only by PatchNeighborBoundary and checks they reproduce a single
# reference grid running the same wall-colliding flow.
#
# NCYCLES is small: stepping splitA then splitB independently each cycle is
# Gauss-Seidel (each sees the other's ghost one cycle staler than vice versa)
# against the reference's Jacobi update, so they only match before that lag
# compounds -- synchronized multi-patch stepping is milestone 6's job.
# ============================================================================
nx2, ny2 = 8, 3   # real (evolving) cells; the reference grid adds its own 1-cell wall-ghost rim on top
dx2, dy2 = 0.1, 0.1
NCYCLES = 3

def makeIC(density, pressure, vx):
    def ic(nodeList, grid):
        rho = nodeList.getFieldDouble("density")
        u = nodeList.getFieldDouble("specificInternalEnergy")
        v = nodeList.getFieldVector2d("velocity")
        for idx in range(grid.size()):
            rho.setValue(idx, density)
            u.setValue(idx, pressure / ((GAMMA - 1) * density))
            v.setValue(idx, Vector2d(vx, 0.0))
    return ic

DENSITY, PRESSURE, VX = 1.0, 1.0, 0.3
setIC = makeIC(DENSITY, PRESSURE, VX)

# Reference: single grid spanning the whole domain. Built (nx2+2) x (ny2+2)
# since a full ReflectingGridBoundary treats the outermost rim as wall-ghost,
# not real cells -- real cells are at local indices 1..nx2 / 1..ny2.
refGrid = Grid2d(nx2 + 2, ny2 + 2, dx2, dy2)
refNodeList = NodeList(refGrid.size())
refPhysics = physicsFactory(refNodeList, refGrid)
refPhysics.addBoundary(ReflectingGridBoundary2d(grid=refGrid))
refIntegrator = RungeKutta2Integrator2d([refPhysics], dtmin=1e-7)
setIC(refNodeList, refGrid)

# --- Split: two patches, joined only by PatchNeighborBoundary. ---
splitA = buildPatch(box=(0, 0, 3, ny2 - 1), level=0, dx=dx2, dy=dy2, nghost=NGHOST,
                     physicsFactory=physicsFactory, integratorClass=RungeKutta2Integrator2d, dtmin=1e-7)
splitB = buildPatch(box=(4, 0, nx2 - 1, ny2 - 1), level=0, dx=dx2, dy=dy2, nghost=NGHOST,
                     physicsFactory=physicsFactory, integratorClass=RungeKutta2Integrator2d, dtmin=1e-7)

leftWall = ReflectingGridBoundary2d(grid=splitA.grid); leftWall.setFaces(["left", "top", "bottom"])
rightWall = ReflectingGridBoundary2d(grid=splitB.grid); rightWall.setFaces(["right", "top", "bottom"])
splitA.physics.addBoundary(leftWall)
splitB.physics.addBoundary(rightWall)
assert wireSiblingBoundaries(splitA, splitB, NGHOST), "expected splitA/splitB to share a face"

setIC(splitA.nodeList, splitA.grid)
setIC(splitB.nodeList, splitB.grid)

# Integrator.Step() doesn't self-initialize; Controller.Step() normally calls Initialize() on cycle 0.
refIntegrator.Initialize()
splitA.integrator.Initialize()
splitB.integrator.Initialize()

for _ in range(NCYCLES):
    refIntegrator.Step()
    splitA.integrator.Step()
    splitB.integrator.Step()

refDensity = refNodeList.getFieldDouble("density")
refPressure = refNodeList.getFieldDouble("pressure")
refVelocity = refNodeList.getFieldVector2d("velocity")

maxErr = 0.0
for j in range(ny2):
    for i in range(nx2):
        refIdx = refGrid.index(i + 1, j + 1, 0)   # +1: skip the reference grid's own wall-ghost rim
        if i < 4:
            patch, li = splitA, i
        else:
            patch, li = splitB, i - 4
        pIdx = patch.grid.index(NGHOST + li, NGHOST + j, 0)

        pDensity = patch.nodeList.getFieldDouble("density")
        pPressure = patch.nodeList.getFieldDouble("pressure")
        pVelocity = patch.nodeList.getFieldVector2d("velocity")

        maxErr = max(maxErr, abs(pDensity[pIdx] - refDensity[refIdx]))
        maxErr = max(maxErr, abs(pPressure[pIdx] - refPressure[refIdx]))
        maxErr = max(maxErr, abs(pVelocity[pIdx].x - refVelocity[refIdx].x))
        maxErr = max(maxErr, abs(pVelocity[pIdx].y - refVelocity[refIdx].y))

print(f"max |split - reference| after {NCYCLES} cycles: {maxErr:.3e}")
assert maxErr < 1e-12, f"split-domain result diverged from single-grid reference: maxErr={maxErr:.3e}"
print("[ok] AMR milestone 2 (PatchNeighborBoundary) passes: split domain reproduces the single-grid reference to near machine precision.")

# ============================================================================
# Milestone 3: AMR::CoarseFineBoundary -- piecewise-constant prolongation from
# a coarse parent into a fine patch's ghost rim. Static 2-level hierarchy: no
# subcycling, no restriction, no reflux yet.
# ============================================================================

coarseNx, coarseNy = 10, 10
cdx, cdy = 0.1, 0.1

coarseGrid = Grid2d(coarseNx, coarseNy, cdx, cdy)
coarseNodes = NodeList(coarseGrid.size())
hier = AMRController(coarseGrid, coarseNodes, physicsFactory,
                      RungeKutta2Integrator2d, dtmin=1e-7,
                      refinementRatio=REFINEMENT_RATIO, nghost=NGHOST)

# A distinct value per coarse cell, small enough not to trip the solver's explode guard.
coarseDensity = coarseNodes.getFieldDouble("density")
coarseEnergy = coarseNodes.getFieldDouble("specificInternalEnergy")
coarseVel = coarseNodes.getFieldVector2d("velocity")
for idx in range(coarseGrid.size()):
    coarseDensity.setValue(idx, 1.0 + 0.01 * idx)
    coarseEnergy.setValue(idx, 1.0)
    coarseVel.setValue(idx, Vector2d(0.0, 0.0))

def containingCoarseCell(pos):
    """Coarse cell whose extent contains `pos`, found by scanning coarse cell centers."""
    for cIdx in range(coarseGrid.size()):
        c = coarseGrid.position(cIdx)
        if abs(c.x - pos.x) <= 0.5 * cdx + 1e-12 and abs(c.y - pos.y) <= 0.5 * cdy + 1e-12:
            return cIdx
    return None

def checkGhostMapping(box, label):
    """Wire a fine patch at `box` to the coarse level, step it once, and verify
    every ghost cell got the value of the coarse cell geometrically containing
    it -- checked from world positions, not the wiring's own index math."""
    patch = hier.buildPatch(level=1, box=box)
    paired = wireCoarseFineBoundary(patch, hier.levels[0][0], REFINEMENT_RATIO)

    nx = (box[2] - box[0] + 1) + 2 * NGHOST
    ny = (box[3] - box[1] + 1) + 2 * NGHOST
    expected = nx * ny - (nx - 2 * NGHOST) * (ny - 2 * NGHOST)
    assert paired == expected, f"{label}: paired {paired} cells, expected {expected} ghosts"

    density = patch.nodeList.getFieldDouble("density")
    energy = patch.nodeList.getFieldDouble("specificInternalEnergy")
    vel = patch.nodeList.getFieldVector2d("velocity")
    for idx in range(patch.grid.size()):
        density.setValue(idx, 1.0)
        energy.setValue(idx, 1.0)
        vel.setValue(idx, Vector2d(0.0, 0.0))

    patch.integrator.Initialize()
    patch.integrator.Step()

    checked = 0
    for lj in range(ny):
        for li in range(nx):
            if (NGHOST <= li < nx - NGHOST) and (NGHOST <= lj < ny - NGHOST):
                continue
            gIdx = patch.grid.index(li, lj, 0)
            cIdx = containingCoarseCell(patch.grid.position(gIdx))
            assert cIdx is not None, f"{label}: ghost ({li},{lj}) has no containing coarse cell"
            assert abs(density[gIdx] - coarseDensity[cIdx]) < 1e-14, (
                f"{label}: ghost ({li},{lj}) got {density[gIdx]}, expected coarse cell {cIdx} value {coarseDensity[cIdx]}")
            checked += 1
    print(f"[ok] {label}: all {checked} ghost cells received the value of the coarse cell containing them.")
    return patch

# Aligned: patch edges fall on coarse cell boundaries. Misaligned: they don't,
# exercising the floor division that maps a fine index onto its coarse cell.
fine = checkGhostMapping((3 * REFINEMENT_RATIO, 3 * REFINEMENT_RATIO,
                          6 * REFINEMENT_RATIO + REFINEMENT_RATIO - 1,
                          6 * REFINEMENT_RATIO + REFINEMENT_RATIO - 1), "aligned patch")
checkGhostMapping((5, 5, 12, 12), "misaligned patch")

fineNx = (fine.box[2] - fine.box[0] + 1) + 2 * NGHOST
fineNy = (fine.box[3] - fine.box[1] + 1) + 2 * NGHOST
fineDensity = fine.nodeList.getFieldDouble("density")
fineEnergy = fine.nodeList.getFieldDouble("specificInternalEnergy")
fineVel = fine.nodeList.getFieldVector2d("velocity")

# --- A uniform state must stay exactly uniform: prolongation into the rim
# must not inject spurious structure at the coarse-fine interface. ---
UNIFORM_RHO, UNIFORM_U, UNIFORM_VX = 1.0, 1.0, 0.3
for idx in range(coarseGrid.size()):
    coarseDensity.setValue(idx, UNIFORM_RHO)
    coarseEnergy.setValue(idx, UNIFORM_U)
    coarseVel.setValue(idx, Vector2d(UNIFORM_VX, 0.0))
for idx in range(fine.grid.size()):
    fineDensity.setValue(idx, UNIFORM_RHO)
    fineEnergy.setValue(idx, UNIFORM_U)
    fineVel.setValue(idx, Vector2d(UNIFORM_VX, 0.0))

for _ in range(10):
    fine.integrator.Step()

worstUniform = 0.0
for lj in range(NGHOST, fineNy - NGHOST):
    for li in range(NGHOST, fineNx - NGHOST):
        idx = fine.grid.index(li, lj, 0)
        worstUniform = max(worstUniform, abs(fineDensity[idx] - UNIFORM_RHO))
        worstUniform = max(worstUniform, abs(fineEnergy[idx] - UNIFORM_U))
        worstUniform = max(worstUniform, abs(fineVel[idx].x - UNIFORM_VX))
        worstUniform = max(worstUniform, abs(fineVel[idx].y))

print(f"max deviation from uniform on the fine patch after 10 cycles: {worstUniform:.3e}")
assert worstUniform < 1e-12, f"coarse-fine interface injected structure into a uniform state: {worstUniform:.3e}"
print("[ok] AMR milestone 3 (CoarseFineBoundary) passes: prolongation maps ghosts correctly and preserves a uniform state.")

# ============================================================================
# Milestone 4: AMR::RestrictionOperator -- fine->coarse averaging. Mass,
# momentum and total energy over the covered region must be identical before
# and after, which a naive arithmetic mean of velocity/sie would not give.
# ============================================================================

rGrid = Grid2d(coarseNx, coarseNy, cdx, cdy)
rNodes = NodeList(rGrid.size())
rHier = AMRController(rGrid, rNodes, physicsFactory, RungeKutta2Integrator2d,
                       dtmin=1e-7, refinementRatio=REFINEMENT_RATIO, nghost=NGHOST)
rCoarse = rHier.levels[0][0]

COVERED_LO, COVERED_HI = 3, 6   # coarse cells covered by the fine patch, aligned so each is covered in full
rFine = rHier.buildPatch(level=1, box=(COVERED_LO * REFINEMENT_RATIO, COVERED_LO * REFINEMENT_RATIO,
                                        COVERED_HI * REFINEMENT_RATIO + REFINEMENT_RATIO - 1,
                                        COVERED_HI * REFINEMENT_RATIO + REFINEMENT_RATIO - 1))
restriction = buildRestriction(rFine, rCoarse, REFINEMENT_RATIO, eos)

# Coarse state starts uniform, with a sentinel outside the covered region to
# confirm restriction leaves uncovered cells alone.
rcRho = rNodes.getFieldDouble("density")
rcU = rNodes.getFieldDouble("specificInternalEnergy")
rcVel = rNodes.getFieldVector2d("velocity")
for idx in range(rGrid.size()):
    rcRho.setValue(idx, 1.0)
    rcU.setValue(idx, 1.0)
    rcVel.setValue(idx, Vector2d(0.0, 0.0))
sentinelIdx = rGrid.index(0, 0, 0)
rcRho.setValue(sentinelIdx, 7.5)

# Non-uniform fine state, so a wrong (unweighted) average would show up.
rfRho = rFine.nodeList.getFieldDouble("density")
rfU = rFine.nodeList.getFieldDouble("specificInternalEnergy")
rfVel = rFine.nodeList.getFieldVector2d("velocity")
fineNx4 = (rFine.box[2] - rFine.box[0] + 1) + 2 * NGHOST
fineNy4 = (rFine.box[3] - rFine.box[1] + 1) + 2 * NGHOST
for lj in range(fineNy4):
    for li in range(fineNx4):
        idx = rFine.grid.index(li, lj, 0)
        rfRho.setValue(idx, 1.0 + 0.3 * ((li * 7 + lj * 3) % 5))
        rfU.setValue(idx, 1.0 + 0.1 * ((li + 2 * lj) % 4))
        rfVel.setValue(idx, Vector2d(0.2 * ((li % 3) - 1), 0.15 * ((lj % 3) - 1)))

def totals(nodeList, grid, cells):
    """(mass, momentum_x, momentum_y, total energy) summed over `cells`."""
    rho = nodeList.getFieldDouble("density")
    u = nodeList.getFieldDouble("specificInternalEnergy")
    vel = nodeList.getFieldVector2d("velocity")
    m = px = py = e = 0.0
    for idx in cells:
        vol = grid.cellVolume(idx)
        dm = rho[idx] * vol
        v = vel[idx]
        m += dm
        px += dm * v.x
        py += dm * v.y
        e += dm * (u[idx] + 0.5 * (v.x * v.x + v.y * v.y))
    return m, px, py, e

fineInterior = [rFine.grid.index(li, lj, 0)
                for lj in range(NGHOST, fineNy4 - NGHOST)
                for li in range(NGHOST, fineNx4 - NGHOST)]
coveredCoarse = [rGrid.index(ci, cj, 0)
                 for cj in range(COVERED_LO, COVERED_HI + 1)
                 for ci in range(COVERED_LO, COVERED_HI + 1)]

fineTotals = totals(rFine.nodeList, rFine.grid, fineInterior)
restriction.apply()
coarseTotals = totals(rNodes, rGrid, coveredCoarse)

labels = ("mass", "momentum x", "momentum y", "total energy")
for label, f, c in zip(labels, fineTotals, coarseTotals):
    scale = max(abs(f), 1e-30)
    relErr = abs(c - f) / scale
    print(f"  {label}: fine={f:.12e}  coarse={c:.12e}  rel err={relErr:.3e}")
    assert relErr < 1e-12, f"restriction did not conserve {label}: rel err {relErr:.3e}"

assert rcRho[sentinelIdx] == 7.5, "restriction modified a coarse cell outside the covered region"
print("[ok] AMR milestone 4 (RestrictionOperator) passes: mass, momentum and energy conserved exactly; uncovered cells untouched.")

# ============================================================================
# Milestone 5: AMR::FluxRegister + Refluxer. Without reflux the coarse and fine
# sides of the interface disagree about how much material crossed it and the
# composite hierarchy's mass drifts; with reflux the interface stops leaking.
# Still a single dt for both levels -- subcycling is milestone 6.
#
# The bar is set by measuring a plain single-grid run of the same problem
# rather than by asserting machine precision: GridHydroBase + a reflecting wall
# is itself only approximately mass-conserving (ApplyBoundaries runs once per
# step, so the wall's ghost cells are stale during RK2's second stage), and
# that pre-existing error is far larger than anything the interface contributes.
# ============================================================================

fxNx, fxNy = 16, 16
fxDx = 0.05
FX_CYCLES = 40
FX_LO, FX_HI = 6, 9

def setBlob(nodeList, grid):
    """Ambient gas with an over-pressured blob straddling the refined region's
    edge, so real flow crosses the coarse-fine interface in both directions."""
    rho = nodeList.getFieldDouble("density")
    u = nodeList.getFieldDouble("specificInternalEnergy")
    vel = nodeList.getFieldVector2d("velocity")
    centre = 0.5 * fxNx * fxDx
    for idx in range(grid.size()):
        pos = grid.position(idx)
        inBlob = (pos.x - centre) ** 2 + (pos.y - centre) ** 2 < (2.5 * fxDx) ** 2
        rho.setValue(idx, 2.0 if inBlob else 1.0)
        u.setValue(idx, (5.0 if inBlob else 1.0) / ((GAMMA - 1) * (2.0 if inBlob else 1.0)))
        vel.setValue(idx, Vector2d(0.0, 0.0))

def runHierarchy(useReflux):
    """Advance a 2-level hierarchy for FX_CYCLES and return its relative mass drift."""
    grid = Grid2d(fxNx, fxNy, fxDx, fxDx)
    nodes = NodeList(grid.size())
    hier = AMRController(grid, nodes, physicsFactory, RungeKutta2Integrator2d,
                          dtmin=1e-8, refinementRatio=REFINEMENT_RATIO, nghost=NGHOST,
                          boundaries=[ReflectingGridBoundary2d(grid=grid)])
    coarse = hier.levels[0][0]
    fine = hier.buildPatch(level=1, box=(FX_LO * REFINEMENT_RATIO, FX_LO * REFINEMENT_RATIO,
                                          FX_HI * REFINEMENT_RATIO + REFINEMENT_RATIO - 1,
                                          FX_HI * REFINEMENT_RATIO + REFINEMENT_RATIO - 1))
    wireCoarseFineBoundary(fine, coarse, REFINEMENT_RATIO)
    restriction = buildRestriction(fine, coarse, REFINEMENT_RATIO, eos)
    reflux = buildReflux(fine, coarse, REFINEMENT_RATIO, eos) if useReflux else None

    setBlob(nodes, grid)
    setBlob(fine.nodeList, fine.grid)

    fnx = (fine.box[2] - fine.box[0] + 1) + 2 * NGHOST
    fny = (fine.box[3] - fine.box[1] + 1) + 2 * NGHOST
    fineInterior = [fine.grid.index(li, lj, 0)
                    for lj in range(NGHOST, fny - NGHOST) for li in range(NGHOST, fnx - NGHOST)]
    # Coarse cells that are real (not the domain's wall-ghost rim) and not covered by the fine patch.
    uncovered = [grid.index(ci, cj, 0)
                 for cj in range(1, fxNy - 1) for ci in range(1, fxNx - 1)
                 if not (FX_LO <= ci <= FX_HI and FX_LO <= cj <= FX_HI)]

    cRho = nodes.getFieldDouble("density")
    fRho = fine.nodeList.getFieldDouble("density")
    def mass():
        return (sum(cRho[i] * grid.cellVolume(i) for i in uncovered)
                + sum(fRho[i] * fine.grid.cellVolume(i) for i in fineInterior))

    coarse.integrator.Initialize()
    fine.integrator.Initialize()
    m0 = mass()
    for cyc in range(FX_CYCLES):
        dt = min(coarse.integrator.dt, fine.integrator.dt)
        # Both levels must advance by the same dt for a single-dt reflux to be
        # consistent; restoreState is the only bound way to set an integrator's dt.
        coarse.integrator.restoreState(cyc, cyc * dt, dt)
        fine.integrator.restoreState(cyc, cyc * dt, dt)
        coarse.integrator.Step()
        fine.integrator.Step()
        if reflux:
            reflux.apply(dt)
        restriction.apply()
    return abs(mass() - m0) / m0

def runSingleGrid():
    """Same problem on a plain single grid -- the scheme's own conservation floor."""
    grid = Grid2d(fxNx, fxNy, fxDx, fxDx)
    nodes = NodeList(grid.size())
    physics = physicsFactory(nodes, grid)
    physics.addBoundary(ReflectingGridBoundary2d(grid=grid))
    integrator = RungeKutta2Integrator2d([physics], dtmin=1e-8)
    setBlob(nodes, grid)
    real = [grid.index(i, j, 0) for j in range(1, fxNy - 1) for i in range(1, fxNx - 1)]
    rho = nodes.getFieldDouble("density")
    mass = lambda: sum(rho[i] * grid.cellVolume(i) for i in real)
    integrator.Initialize()
    m0 = mass()
    for _ in range(FX_CYCLES):
        integrator.Step()
    return abs(mass() - m0) / m0

driftBaseline = runSingleGrid()
driftReflux = runHierarchy(useReflux=True)
driftNoReflux = runHierarchy(useReflux=False)

print(f"mass drift over {FX_CYCLES} cycles:  single grid={driftBaseline:.3e}"
      f"  hierarchy+reflux={driftReflux:.3e}  hierarchy no reflux={driftNoReflux:.3e}")
assert driftReflux < max(driftBaseline, 1e-12), (
    f"reflux hierarchy ({driftReflux:.3e}) drifts more than a plain single grid ({driftBaseline:.3e})")
assert driftNoReflux > 10 * driftReflux, (
    f"reflux made no appreciable difference: {driftNoReflux:.3e} vs {driftReflux:.3e}")
print("[ok] AMR milestone 5 (FluxRegister/reflux) passes: the interface stops leaking mass, "
      "leaving the hierarchy no worse than a single grid.")
