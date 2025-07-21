from yggdrasil import *
import numpy as np
import random
from Physics import MagField3d
from Animation import AnimateScatter

class Projected2DView:
    def __init__(self, vector3d_list, axes=(0, 1)):
        self.data = vector3d_list
        self.axis_x, self.axis_y = axes

    def size(self):
        return self.data.size()

    def __getitem__(self, i):
        v = self.data[i]
        class View2D:
            def __init__(self, v, axx, axy):
                self.x = v.x
                self.y = v.y
        return View2D(v, self.axis_x, self.axis_y)

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                       numNodes = 1,
                                       b = 1.0)


    constants  = CGS()
    myNodeList = NodeList(numNodes)
    bVec       = Vector3d(0, 0, b)

    mag = MagField3d(myNodeList, constants, bVec)
    packages = [mag]

    integrator = RungeKutta2Integrator3d(packages=packages, dtmin=1e-5,verbose=False)
    print(integrator)

    mass = myNodeList.mass
    for i in range(numNodes):
        mass[i] = 0.2

    print("numNodes =", myNodeList.numNodes)
    print("field names =", myNodeList.fieldNames)

    myNodeList.position[0]  = Vector3d(1, 0, 0)
    myNodeList.velocity[0]  = Vector3d(0, -5, 0)
    myNodeList.charge[0]    = 1.0

    pos = myNodeList.position
    projected_positions = Projected2DView(pos, axes=(0, 1))  # XY plane

    controller = Controller(integrator=integrator, periodicWork=[], statStep=1)

    bounds = (-10, 10, -10, 10)
    AnimateScatter(bounds, stepper=controller, positions=projected_positions, frames=100, interval=50)
