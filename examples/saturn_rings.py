import numpy as np
import random
from yggdrasil import *
from Animation import *
from math import sin,cos,atan2
from Physics import PointSourceGravity2d, Kinematics2d


if __name__ == "__main__":
    commandLine = CommandLineArguments(nBodies = 200,
                                       cmass = 1.0,
                                       cycles = 1200000)
    
    constants = PhysicalConstants(58232e+3,     # saturn radius in m
                                  5.683e+26,    # saturn mass in kg 
                                  3.156e7       # 1 year in s
                                  ) 
    print("G: ",constants.G)

    cmLoc = Vector2d(0,0)
    mmLoc = Vector2d(60000,-500000)
    mmVel = Vector2d(0,0.1)


    
    from PoissonNodeGenerator import PoissonDisk2d
    pd = PoissonDisk2d(nBodies)

    for i in range(len(pd.positions)):
        pd.positions[i] = (pd.positions[i][0]*6.87,pd.positions[i][1]*6.87)
    pd.positions = [position for position in pd.positions if np.sqrt(position[0]**2 + position[1]**2) > 1.11]

    nBodies = len(pd.positions)
    myNodeList = NodeList(nBodies)

    kinematics = Kinematics2d(nodeList=myNodeList, constants=constants)
    sourceGrav = PointSourceGravity2d(nodeList=myNodeList,
                                      constants=constants,
                                      pointSourceLocation=cmLoc,
                                      pointSourceMass=cmass,
                                      pointSourceVelocity = Vector2d(0,0))
    
    packages = [kinematics, sourceGrav]
    integrator = RungeKutta4Integrator2d(packages=packages,
                                         dtmin=1e-5,verbose=False)


    pos = myNodeList.position
    velocity = myNodeList.velocity

    for i in range(nBodies):
        # r = float(random.uniform(1.11,6.87))
        # t = random.random()*2.0*np.pi
        # pos[i].x = r*cos(t)
        # pos[i].y = r*sin(t)
        pos[i].x = pd.positions[i][0]
        pos[i].y = pd.positions[i][1]
        t = atan2(pos[i].y,pos[i].x)
        r = np.sqrt(pos[i].x**2 + pos[i].y**2)
        v = np.sqrt(cmass*constants.G/r) # 80% of circular velocity
        velocity[i].x = -v*sin(t)
        velocity[i].y = v*cos(t)

    periodicWork = []


    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=10000)

    bounds = (-10, 10, -10, 10)
    AnimateScatter(bounds, stepper=controller, positions=pos, frames=10, interval=50)