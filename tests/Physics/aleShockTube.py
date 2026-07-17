import os
from yggdrasil import *
from Mesh import ALEMesh2d
from Physics import ALEMeshHydroHLLE2d, ALEMeshHydroHLLC2d
from EOS import IdealGasEOS
from Boundaries import ReflectingALEMeshBoundary2d, DirichletALEMeshBoundary2d, OutflowALEMeshBoundary2d
from Utilities import SiloDump
import math

# Sod-shock-tube-style test on a hand-built unstructured strip (not a Grid):
# a thin row of quad cells, Dirichlet (fixed high-pressure reservoir) on the
# left end, Outflow on the right end, Reflecting on the two long channel
# walls. Exercises ALEMeshHydro end to end: face construction (interior +
# boundary), all three ALEMeshBoundary types, and the flux*area/volume
# divergence in ALEMeshHydroBase.
#
# HLLE (two-wave, no exact contact resolution) runs cleanly for hundreds of
# cycles on this problem. HLLC (exact contact resolution, but no sonic
# entropy fix) is limited to a short stable window before failing at the
# leading edge of the rarefaction fan -- a known HLL-family weakness at sonic
# points, not specific to this mesh. Swap the two lines below to see it.

if __name__ == "__main__":
    commandLine = CommandLineArguments(nx = 40,
                                        ny = 2,
                                        dx = 1.0,
                                        dy = 1.0,
                                        cycles = 400,
                                        dumpCycle = 20,
                                        siloDump = False,
                                        rootName = "aleShockTube",
                                        vizDir = "viz")

    if siloDump:
        vizDir = os.path.join(rootName, vizDir)
        os.makedirs(vizDir, exist_ok=True)

    mesh = ALEMesh2d()

    def nodeId(i, j):
        return i + j * (nx + 1)

    def cellId(i, j):
        return i + j * nx

    for j in range(ny + 1):
        for i in range(nx + 1):
            mesh.addNode(Vector2d(float(i) * dx, float(j) * dy))

    for j in range(ny):
        for i in range(nx):
            mesh.addCell([nodeId(i, j), nodeId(i + 1, j), nodeId(i + 1, j + 1), nodeId(i, j + 1)])

    mesh.computeFaces()

    constants = MKS()
    eos = IdealGasEOS(1.4, constants)

    numCells = nx * ny
    myNodeList = NodeList(numCells)
    hydro = ALEMeshHydroHLLE2d(myNodeList, constants, eos, mesh)
    # hydro = ALEMeshHydroHLLC2d(myNodeList, constants, eos, mesh)

    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")

    rhoL, uL_ = 1.0, 2.5
    rhoR, uR_ = 0.125, 2.0

    for j in range(ny):
        for i in range(nx):
            idx = cellId(i, j)
            if i < nx // 2:
                density.setValue(idx, rhoL)
                energy.setValue(idx, uL_)
            else:
                density.setValue(idx, rhoR)
                energy.setValue(idx, uR_)
            velocity.setValue(idx, Vector2d(0.0, 0.0))

    # Classify boundary faces by centroid position (left end / right end /
    # top+bottom channel walls) so each ALEMeshBoundary governs the right set.
    faces = mesh.getFaces()
    leftFaces, rightFaces, wallFaces = [], [], []
    xmax = nx * dx
    for idx, f in enumerate(faces):
        if not f.isBoundary:
            continue
        cx = f.centroid.x
        if abs(cx - 0.0) < 1e-9:
            leftFaces.append(idx)
        elif abs(cx - xmax) < 1e-9:
            rightFaces.append(idx)
        else:
            wallFaces.append(idx)

    print(f"boundary faces: left={len(leftFaces)} right={len(rightFaces)} walls={len(wallFaces)}")
    assert len(leftFaces) == ny
    assert len(rightFaces) == ny
    assert len(wallFaces) == 2 * nx

    inflow = DirichletALEMeshBoundary2d(mesh, rhoL, Vector2d(0.0, 0.0), uL_)
    inflow.setFaces(leftFaces)
    hydro.addBoundary(inflow)

    outflow = OutflowALEMeshBoundary2d(mesh)
    outflow.setFaces(rightFaces)
    hydro.addBoundary(outflow)

    walls = ReflectingALEMeshBoundary2d(mesh)
    walls.setFaces(wallFaces)
    hydro.addBoundary(walls)

    integrator = RungeKutta4Integrator2d([hydro], dtmin=1e-6, verbose=False)

    periodicWork = []
    if siloDump:
        # mesh= (not grid=) writes the ALEMesh's real polygon topology --
        # arbitrary cell shapes, not just quads -- rather than a point cloud
        # at cell centroids. Open the resulting cycle=NNN.silo files in VisIt.
        dump = SiloDump(baseName=os.path.join(vizDir, rootName),
                         nodeList=myNodeList,
                         fieldNames=["density", "specificInternalEnergy", "pressure", "velocity"],
                         dumpCycle=dumpCycle,
                         mesh=mesh)
        periodicWork += [dump]

    controller = Controller(integrator=integrator, periodicWork=periodicWork, statStep=50)
    controller.Step(cycles)

    profile = [round(density[cellId(i, 0)], 3) for i in range(nx)]
    print("density profile along the strip:")
    print(profile)

    bad = any(not math.isfinite(density[idx]) for idx in range(numCells))
    assert not bad, "non-finite density values"

    # A real shock/rarefaction structure should have developed: density
    # should stay roughly within the initial [rhoR, rhoL] bounds (no wild
    # overshoot/blowup) and should no longer be the flat two-block step it
    # started as.
    assert max(profile) <= rhoL + 0.05, f"density overshoot: max={max(profile)}"
    assert min(profile) >= rhoR - 0.05, f"density undershoot: min={min(profile)}"
    distinctValues = len(set(profile))
    print("distinct density values along the strip:", distinctValues)
    assert distinctValues > 5, "expected a developed shock/rarefaction structure, not a flat step"

    print("PASS: aleShockTube")
