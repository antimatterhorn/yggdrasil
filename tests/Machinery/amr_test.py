from yggdrasil import *
from Mesh import Grid2d
from EOS import IdealGasEOS
from Physics import GridHydroHLLC2d
from Boundaries import ReflectingGridBoundary2d
from AMRController import AMRController, buildPatch, wireSiblingBoundaries

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
