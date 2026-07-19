import os
from yggdrasil import *
from Mesh import ALEMesh2d
from Physics import ALEMeshHydroHLLE2d
from EOS import IdealGasEOS
from Boundaries import DirichletALEMeshBoundary2d
from Utilities import SiloDump
import math

# Mesh-motion + remap sanity check: a uniform flow (same density/pressure/
# velocity everywhere) on a mesh whose every node is prescribed to move at
# exactly the flow's own velocity -- the whole mesh rigidly translates along
# with the material. Physically nothing is happening: in the frame co-moving
# with the mesh the fluid is perfectly static, so density/pressure/velocity
# must stay exactly (to floating-point roundoff) unchanged for the entire
# run, no matter how far the mesh visibly translates. This exercises the ALE
# flux (S = w.n != 0 at every boundary face along the translation direction),
# the FinalizeStep motion hook, and the extensive-quantity remap together --
# and is exactly the kind of case where a subtle ALE-flux term/sign error
# (see HLL.cc's derivation comments -- the naive "shift vn everywhere"
# shortcut is wrong specifically for the energy term) would show up as
# density/pressure drift rather than staying silent.
#
# Dirichlet boundaries on all four sides, matching the uniform interior state
# exactly, sidestep a separate, currently-unaddressed question (out of scope
# for this pass): whether ReflectingALEMeshBoundary's ghostState needs to be
# face-velocity-aware to correctly handle a wall that is itself moving.
# Dirichlet's ghost state doesn't depend on face motion at all, so it still
# genuinely exercises the ALE flux without depending on that answer.

if __name__ == "__main__":
    commandLine = CommandLineArguments(nx=10,
                                        ny=2,
                                        dx=1.0,
                                        dy=1.0,
                                        cycles=250,
                                        vx0=0.3,
                                        dumpCycle=20,
                                        siloDump=False,
                                        rootName="aleRigidTranslation",
                                        vizDir="viz")

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
    gamma = 1.4
    eos = IdealGasEOS(gamma, constants)

    numCells = nx * ny
    myNodeList = NodeList(numCells)

    rho0 = 1.0
    u0 = 2.0
    p0 = (gamma - 1.0) * rho0 * u0
    v0 = Vector2d(vx0, 0.0)

    # Node-centered mesh-velocity NodeList: every node moves at exactly v0.
    numNodes = len(mesh.getNodes())
    nodeVelocities = NodeList(numNodes)

    # ALEMeshHydroHLLE2d's constructor enrolls density/velocity/
    # specificInternalEnergy/pressure/soundSpeed on myNodeList and "velocity"
    # on nodeVelocities -- construct it before fetching any of those fields.
    hydro = ALEMeshHydroHLLE2d(myNodeList, constants, eos, mesh, nodeVelocities)

    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")
    for j in range(ny):
        for i in range(nx):
            idx = cellId(i, j)
            density.setValue(idx, rho0)
            energy.setValue(idx, u0)
            velocity.setValue(idx, v0)

    w = nodeVelocities.getFieldVector2d("velocity")
    for n in range(numNodes):
        w.setValue(n, v0)

    faces = mesh.getFaces()
    boundaryFaceIds = [idx for idx, f in enumerate(faces) if f.isBoundary]
    print("boundary faces:", len(boundaryFaceIds))
    assert len(boundaryFaceIds) == 2 * nx + 2 * ny

    allBoundary = DirichletALEMeshBoundary2d(mesh, rho0, v0, u0)
    allBoundary.setFaces(boundaryFaceIds)
    hydro.addBoundary(allBoundary)

    integrator = RungeKutta4Integrator2d([hydro], dtmin=1e-6, verbose=False)

    periodicWork = []
    if siloDump:
        # mesh= (not grid=) writes the ALEMesh's real polygon topology, so
        # the mesh's visible translation shows up across successive dumps.
        dump = SiloDump(baseName=os.path.join(vizDir, rootName),
                         nodeList=myNodeList,
                         fieldNames=["density", "specificInternalEnergy", "pressure", "velocity"],
                         dumpCycle=dumpCycle,
                         mesh=mesh)
        periodicWork += [dump]

    controller = Controller(integrator=integrator, periodicWork=periodicWork, statStep=25)

    initialNodePositions = [Vector2d(p.x, p.y) for p in mesh.getNodes()]

    controller.Step(cycles)

    elapsed = controller.time
    print("elapsed time:", elapsed)

    # Density/pressure/velocity should be unchanged (to roundoff) everywhere.
    pressure = myNodeList.getFieldDouble("pressure")
    maxRhoDrift = 0.0
    maxPressureDrift = 0.0
    maxVelDrift = 0.0
    for idx in range(numCells):
        maxRhoDrift = max(maxRhoDrift, abs(density[idx] - rho0) / rho0)
        maxPressureDrift = max(maxPressureDrift, abs(pressure[idx] - p0) / p0)
        vi = velocity[idx]
        maxVelDrift = max(maxVelDrift, abs(vi.x - v0.x), abs(vi.y - v0.y))

    print("max relative density drift:", maxRhoDrift)
    print("max relative pressure drift:", maxPressureDrift)
    print("max velocity component drift:", maxVelDrift)

    assert math.isfinite(maxRhoDrift) and maxRhoDrift < 1e-6, f"density drifted: {maxRhoDrift}"
    assert math.isfinite(maxPressureDrift) and maxPressureDrift < 1e-6, f"pressure drifted: {maxPressureDrift}"
    assert math.isfinite(maxVelDrift) and maxVelDrift < 1e-6, f"velocity drifted: {maxVelDrift}"

    # The mesh must have genuinely moved -- this isn't testing "nothing
    # happened" by accident. Every node should have translated by v0*elapsed.
    maxPosError = 0.0
    finalNodePositions = mesh.getNodes()
    for n in range(numNodes):
        expected = initialNodePositions[n] + v0 * elapsed
        actual = finalNodePositions[n]
        maxPosError = max(maxPosError, abs(actual.x - expected.x), abs(actual.y - expected.y))
    print("max node position error vs v0*elapsed:", maxPosError)

    expectedTravel = vx0 * elapsed
    print("expected travel distance:", expectedTravel)
    assert expectedTravel > dx, f"mesh barely moved ({expectedTravel}) -- not a meaningful motion test"
    assert maxPosError < 1e-9, f"node positions don't match v0*elapsed: {maxPosError}"

    print("PASS: aleRigidTranslation")
