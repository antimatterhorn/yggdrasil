from yggdrasil import *
from Physics import NBodyGravity2d, TreeGravity2d, FMMGravity2d, Kinematics2d
import numpy as np
import math

def totalEnergy(numNodes, pos, vel, masses, G, eps2):
    ke = sum(0.5 * masses[i] * (vel[i][0]**2 + vel[i][1]**2) for i in range(numNodes))
    pe = 0.0
    for i in range(numNodes):
        for j in range(i + 1, numNodes):
            dx = pos[i][0] - pos[j][0]
            dy = pos[i][1] - pos[j][1]
            pe -= G * masses[i] * masses[j] / math.sqrt(dx*dx + dy*dy + eps2)
    return ke + pe


def runTrajectory(GravityClass, gravityKwargs, numNodes, pos0, vel0, masses, constants,
                   dtmin, nsteps, sampleEvery, eps2):
    nodeList = NodeList(numNodes)
    kinematics = Kinematics2d(nodeList=nodeList, constants=constants)
    gravity = GravityClass(nodeList=nodeList, constants=constants, **gravityKwargs)

    posField = nodeList.getFieldVector2d("position")
    velField = nodeList.getFieldVector2d("velocity")
    massField = nodeList.getFieldDouble("mass")
    for i in range(numNodes):
        posField.setValue(i, Vector2d(*pos0[i]))
        velField.setValue(i, Vector2d(*vel0[i]))
        massField.setValue(i, masses[i])

    integrator = RungeKutta4Integrator2d(packages=[kinematics, gravity], dtmin=dtmin, verbose=False)
    controller = Controller(integrator=integrator, periodicWork=[], statStep=10**9, tstop=1e30)

    energies = []
    for _ in range(nsteps // sampleEvery):
        controller.Step(nsteps=sampleEvery)
        pos = [(posField[i].x, posField[i].y) for i in range(numNodes)]
        vel = [(velField[i].x, velField[i].y) for i in range(numNodes)]
        energies.append(totalEnergy(numNodes, pos, vel, masses, constants.G, eps2))
    # nodeList must outlive posField/velField -- kept alive by this function's own scope until return.
    return energies


if __name__ == "__main__":
    GMKS = 6.67430e-11
    constants = PhysicalConstants(1.0, 1.0, 1.0 / math.sqrt(GMKS), 1.0, 1.0)  # G == 1 in these units

    rng = np.random.default_rng(11)
    numNodes = 40
    minSep = 1.0
    pos0 = []
    while len(pos0) < numNodes:
        candidate = rng.uniform(-8.0, 8.0, size=2)
        if all(np.linalg.norm(candidate - p) >= minSep for p in pos0):
            pos0.append(candidate)
    pos0 = np.array(pos0)
    vel0 = rng.uniform(-0.01, 0.01, size=(numNodes, 2))
    masses = rng.uniform(0.5, 2.0, size=numNodes)
    eps2 = 0.25
    dtmin, nsteps, sampleEvery = 0.02, 400, 20

    eNBody = runTrajectory(NBodyGravity2d, dict(plummerLength=eps2),
        numNodes, pos0, vel0, masses, constants, dtmin, nsteps, sampleEvery, eps2)
    eTree = runTrajectory(TreeGravity2d, dict(plummerLength=eps2),
        numNodes, pos0, vel0, masses, constants, dtmin, nsteps, sampleEvery, eps2)
    eFMM = runTrajectory(FMMGravity2d, dict(plummerLength=eps2, maxSourcesPerLeaf=16),
        numNodes, pos0, vel0, masses, constants, dtmin, nsteps, sampleEvery, eps2)

    # Chaotic close encounters dominate absolute drift for every solver here, so the real check is
    # whether FMMGravity's trajectory tracks NBodyGravity's, not absolute energy conservation.
    relDiff = [abs(eFMM[k] - eNBody[k]) / abs(eNBody[k]) for k in range(len(eNBody))]
    maxRelDiff = max(relDiff)

    # Context only, not asserted on: TreeGravity's monopole force is expected to diverge from the
    # exact trajectory more than FMMGravity's quadrupole-corrected one, via chaotic amplification.
    relDiffTree = [abs(eTree[k] - eNBody[k]) / abs(eNBody[k]) for k in range(len(eNBody))]

    print(f"FMMGravity2d vs NBodyGravity2d: max relative energy difference over trajectory = {maxRelDiff:.3e}")
    print(f"TreeGravity2d vs NBodyGravity2d (context only): max relative energy difference = {max(relDiffTree):.3e}")

    assert maxRelDiff < 0.10, f"FMMGravity2d's trajectory diverged from NBodyGravity2d's: {maxRelDiff:.3e}"
    print("PASS")
