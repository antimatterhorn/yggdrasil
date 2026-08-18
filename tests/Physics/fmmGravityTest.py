from yggdrasil import *
from Physics import (NBodyGravity1d, FMMGravity1d, Kinematics1d,
                      NBodyGravity2d, FMMGravity2d, Kinematics2d,
                      NBodyGravity3d, FMMGravity3d, Kinematics3d)
import Integrators
import numpy as np

def buildNodeList(dim, numNodes, positions, masses, constants, Kinematics, gravityFactory):
    nodeList = NodeList(numNodes)
    kinematics = Kinematics(nodeList=nodeList, constants=constants)
    gravity = gravityFactory(nodeList, constants)

    VectorCtor = globals()[f"Vector{dim}d"]
    posField = getattr(nodeList, f"getFieldVector{dim}d")("position")
    velField = getattr(nodeList, f"getFieldVector{dim}d")("velocity")
    massField = nodeList.getFieldDouble("mass")
    for i in range(numNodes):
        posField.setValue(i, VectorCtor(*positions[i]))
        velField.setValue(i, VectorCtor(*([0.0] * dim)))
        massField.setValue(i, masses[i])

    integrator = getattr(Integrators, f"RungeKutta4Integrator{dim}d")(
        packages=[kinematics, gravity], dtmin=1e-8, verbose=False)
    controller = Controller(integrator=integrator, periodicWork=[], statStep=1, tstop=1e30)
    controller.Step(nsteps=1)
    # nodeList must outlive the returned Field -- it isn't kept alive by the Field itself.
    return nodeList, getattr(nodeList, f"getFieldVector{dim}d")("acceleration")


def checkDim(dim, NBody, FMM, Kinematics, numNodes=300, plummerLength=1e-6, tol=0.05):
    rng = np.random.default_rng(3)
    positions = rng.uniform(-10.0, 10.0, size=(numNodes, dim))
    masses = rng.uniform(0.5, 2.0, size=numNodes)
    constants = PhysicalConstants(1.0, 1.0, 1.0, 1.0, 1.0)

    nlNBody, accNBody = buildNodeList(dim, numNodes, positions, masses, constants, Kinematics,
        lambda nl, c: NBody(nodeList=nl, constants=c, plummerLength=plummerLength))
    nlFMM, accFMM = buildNodeList(dim, numNodes, positions, masses, constants, Kinematics,
        lambda nl, c: FMM(nodeList=nl, constants=c, plummerLength=plummerLength, maxSourcesPerLeaf=16))

    diffSq = 0.0
    exactSq = 0.0
    for i in range(numNodes):
        d = accFMM[i] - accNBody[i]
        diffSq += d.mag2
        exactSq += accNBody[i].mag2

    l2RelErr = (diffSq / exactSq) ** 0.5
    print(f"dim={dim}  L2 relative error (FMMGravity vs NBodyGravity): {l2RelErr:.3e}")
    assert l2RelErr < tol, f"FMMGravity{dim}d diverged too far from NBodyGravity{dim}d: L2 rel err {l2RelErr:.3e}"


if __name__ == "__main__":
    checkDim(1, NBodyGravity1d, FMMGravity1d, Kinematics1d)
    checkDim(2, NBodyGravity2d, FMMGravity2d, Kinematics2d)
    checkDim(3, NBodyGravity3d, FMMGravity3d, Kinematics3d)
    print("PASS")
