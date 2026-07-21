from yggdrasil import *
from Physics import TreeGravity2d, Kinematics2d
from RandomNodeGenerator import RandomNodeGenerator2d
import random

def run():
    numNodes = 20
    cycles = 100

    random.seed(42)

    bounds  = [[-1, -1], [1, 1]]
    vbounds = [[-0.015, -0.015], [0.015, 0.015]]
    posGenerator = RandomNodeGenerator2d(numNodes=numNodes, bounds=bounds)
    velGenerator = RandomNodeGenerator2d(numNodes=numNodes, bounds=vbounds)

    myNodeList = NodeList(numNodes)

    constants = PhysicalConstants(6.378e+6,     # earth radius in m
                                  5.972e+24,    # earth mass in kg
                                  1.0,          # s
                                  1.0,
                                  1.0)

    kinematics = Kinematics2d(nodeList=myNodeList, constants=constants)
    treeGrav = TreeGravity2d(nodeList=myNodeList,
                              constants=constants,
                              plummerLength=0.01)
    packages = [kinematics, treeGrav]

    positions = myNodeList.position
    velocity  = myNodeList.velocity
    mass      = myNodeList.mass
    for i in range(numNodes):
        mass.setValue(i, 10.0)
        positions.setValue(i, Vector2d(posGenerator.positions[i][0], posGenerator.positions[i][1]))
        velocity.setValue(i, Vector2d(velGenerator.positions[i][0], velGenerator.positions[i][1]) * 0.5)

    integrator = RungeKutta4Integrator2d(packages=packages, dtmin=0.01, verbose=False)

    controller = Controller(integrator=integrator, periodicWork=[], statStep=100000)

    controller.Step(cycles)

    out = []
    out.append(integrator.time)
    for i in range(numNodes):
        out.append(positions[i].x)
        out.append(positions[i].y)
    return {"mode": "snapshot", "values": out}

if __name__ == "__main__":
    print(run())
