from yggdrasil import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import random
from Physics import ConstantGravity2d, Kinetics2d
#from Boundaries import SphereCollider2d, BoxCollider2d, MotionConstraint2d
from HCPNodeGenerator import HCPNodeGenerator2d
from LatticeNodeGenerator import Lattice2d
from Animation import AnimateScatter

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                       nx = 10,
                                       ny = 10,
                                       g = -9.8)


    constants = MKS()
    numNodes = nx*ny
    myNodeList = NodeList(numNodes)
    gravVec = Vector2d(0, g)

    lattice = HCPNodeGenerator2d(nx,ny)
    print("hcp gen")


    constantGravity = ConstantGravity2d(myNodeList, constants, gravVec)
    kinetics = Kinetics2d(myNodeList,constants)
    packages = [constantGravity,kinetics]

    for i in range(numNodes):
        lattice.positions[i] = [lattice.positions[i][0]*10,lattice.positions[i][1]+8]

    myNodeList.updatePositions2d(lattice.positions)

    print("nodelist updated")  

    colliders = []
    cbounds = []
    spacing_x = 2  # Horizontal spacing between colliders
    spacing_y = (6 - 2) / 5  # Vertical spacing between colliders to span from 2 to 8
    rows = 8
    collider_radius = 0.3


    # mc = MotionConstraint2d(myNodeList,[0,1],Vector2d(1,1))
    # cbounds.append(mc)


    integrator = RungeKutta2Integrator2d(packages=packages, dtmin=0.01,verbose=False)
    print(integrator)


    rad = myNodeList.getFieldDouble("radius")
    mass = myNodeList.getFieldDouble("mass")
    for i in range(numNodes):
        rad.setValue(i, 0.2)
        mass.setValue(i,0.2)

    print("numNodes =", myNodeList.numNodes)
    print("field names =", myNodeList.fieldNames)

    pos = myNodeList.getFieldVector2d("position")
    print(pos[5])


    controller = Controller(integrator=integrator, periodicWork=[], statStep=1)

    bounds = (-10, 10, 0, 10)
    AnimateScatter(bounds, stepper=controller, positions=pos, frames=100, interval=50)
