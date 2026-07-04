import numpy as np
import os
import random
from yggdrasil import *
from Animation import *
from math import sin,cos,atan2
from Physics import PointSourceGravity2d, Kinematics2d
from IO import RestartWriter


def keplerPosition(time, GM, a, e):
    """Position (x, y) of a body on a fixed Kepler ellipse of semi-major axis
    a and eccentricity e around a mass at the origin, at the given time."""
    n = np.sqrt(GM / a**3)
    M = n * time
    E = M
    for _ in range(50):
        E -= (E - e*sin(E) - M) / (1 - e*cos(E))
    x = a * (cos(E) - e)
    y = a * np.sqrt(1 - e**2) * sin(E)
    return x, y


class MoonMover:
    """Drives a PointSourceGravity's location along a fixed Kepler ellipse
    around Saturn each cycle, so it pulls on the ring particles from the
    moon's true orbital position without simulating the moon's own dynamics."""
    def __init__(self, moonGravity, GM, a, e, workCycle=1):
        self.moonGravity = moonGravity
        self.GM = GM
        self.a = a
        self.e = e
        self.cycle = workCycle

    def __call__(self, cycle, time, dt):
        x, y = keplerPosition(time, self.GM, self.a, self.e)
        self.moonGravity.pointSourceLocation = Vector2d(x, y)


if __name__ == "__main__":
    commandLine = CommandLineArguments(nBodies = 2000,
                                        cmass = 1.0,
                                        mmass = 2e-3,
                                        moonSemiMajorAxis = 9,
                                        moonEccentricity = 0.02,
                                        cycles = 1200000,
                                        rootName = "saturn-rings",
                                        restartDir = "sr/restart",
                                        restartCycle = 500,
                                        restoreCycle = None)
    os.makedirs(restartDir, exist_ok=True)

    constants = PhysicalConstants(58232e+3,     # saturn radius in m
                                  5.683e+26,    # saturn mass in kg
                                  3.156e7       # 1 year in s
                                  )
    print("G: ",constants.G)

    cmLoc = Vector2d(0,0)

    from PoissonNodeGenerator import PoissonDisk2d
    pd = PoissonDisk2d(nBodies,seed=1)

    for i in range(len(pd.positions)):
        pd.positions[i] = (pd.positions[i][0]*6.87,pd.positions[i][1]*6.87)
    pd.positions = [position for position in pd.positions if np.sqrt(position[0]**2 + position[1]**2) > 1.11]

    nBodies = len(pd.positions)
    myNodeList = NodeList(nBodies)

    kinematics = Kinematics2d(nodeList=myNodeList, constants=constants)
    saturnGrav = PointSourceGravity2d(nodeList=myNodeList,
                                      constants=constants,
                                      pointSourceLocation=cmLoc,
                                      pointSourceMass=cmass,
                                      pointSourceVelocity = Vector2d(0,0))

    moonGM = cmass * constants.G
    moonX, moonY = keplerPosition(0.0, moonGM, moonSemiMajorAxis, moonEccentricity)
    moonGrav = PointSourceGravity2d(nodeList=myNodeList,
                                    constants=constants,
                                    pointSourceLocation=Vector2d(moonX, moonY),
                                    pointSourceMass=mmass,
                                    pointSourceVelocity = Vector2d(0,0))

    packages = [kinematics, saturnGrav, moonGrav]
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

    moonMover = MoonMover(moonGrav, moonGM, moonSemiMajorAxis, moonEccentricity, workCycle=1)
    periodicWork = [moonMover]

    restoreIfAvailable(myNodeList, integrator, restartDir, rootName, restoreCycle)
    restartWriter = RestartWriter(myNodeList, integrator)
    def dropRestart(cycle, time, dt):
        restartWriter.write(restartFileName(restartDir, rootName, cycle))
    dropRestart.cycle = restartCycle
    periodicWork += [dropRestart]


    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=10000)

    bounds = (-10, 10, -10, 10)
    AnimateScatter(bounds, stepper=controller, positions=pos, frames=10, interval=50,
                  extra_points=lambda: [(cmLoc.x, cmLoc.y),
                                         (moonGrav.pointSourceLocation.x, moonGrav.pointSourceLocation.y)],
                  extra_colors=['orange', 'red'], extra_size=80)