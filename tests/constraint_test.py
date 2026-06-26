from yggdrasil import *
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import random
from Physics import ConstantGravity2d, Kinetics2d
from Boundaries import SphereCollider2d, BoxCollider2d, MotionConstraint2d
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

    lattice         = HCPNodeGenerator2d(nx,ny)
    constantGravity = ConstantGravity2d(myNodeList, constants, gravVec)
    kinetics        = Kinetics2d(myNodeList,constants)
    packages        = [constantGravity,kinetics]

    for i in range(numNodes):
        lattice.positions[i] = [lattice.positions[i][0]*8,lattice.positions[i][1]+8]

    myNodeList.updatePositions2d(lattice.positions)
    cbounds = []
    indices = []
    pos = myNodeList.getFieldVector2d("position")
    rad = myNodeList.getFieldDouble("radius")
    mass = myNodeList.getFieldDouble("mass")
    for i in range(numNodes):
        rad.setValue(i, 0.2)
        mass.setValue(i,0.2)
    

    for i in range(numNodes):
        if pos[i].x < -7.5:
            indices.append(i)

    mc = MotionConstraint2d(myNodeList,indices,Vector2d(1,1))
    cbounds.append(mc)

    for bound in cbounds:
        constantGravity.addBoundary(bound)

    integrator = RungeKutta2Integrator2d(packages=packages, dtmin=0.01,verbose=False)



    controller = Controller(integrator=integrator, periodicWork=[], statStep=1)

    bounds = (-10, 10, 0, 10)
    AnimateScatter(bounds, stepper=controller, positions=pos, frames=100, interval=50)
