import numpy as np
import os
import random
from yggdrasil import *
from math import sin,cos,atan2
from Physics import PointSourceGravity2d, Kinematics2d
from IO import RestartWriter
from Utilities import SiloDump


def keplerPosition(time, GM, a, e, M0=0.0, retrograde=False):
    """Position (x, y) of a body on a fixed Kepler ellipse of semi-major axis
    a and eccentricity e around a mass at the origin, at the given time.
    M0 sets the orbital phase at time=0; retrograde flips the direction of
    travel."""
    n = np.sqrt(GM / a**3)
    if retrograde:
        n = -n
    M = n * time + M0
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
    def __init__(self, moonGravity, GM, a, e, M0=0.0, retrograde=False, workCycle=1):
        self.moonGravity = moonGravity
        self.GM = GM
        self.a = a
        self.e = e
        self.M0 = M0
        self.retrograde = retrograde
        self.cycle = workCycle

    def __call__(self, cycle, time, dt):
        x, y = keplerPosition(time, self.GM, self.a, self.e, self.M0, self.retrograde)
        self.moonGravity.pointSourceLocation = Vector2d(x, y)


# Orbital elements and mass of Saturn's major inner/mid moons, in SI units
# (meters, kg) -- converted to simulation units from PhysicalConstants below.
# All seven of these orbit prograde; retrograde is kept per-moon for
# generality (Saturn's retrograde moons are irregular outer ones like Phoebe,
# not any of these).
SATURN_MOONS = {
    "Mimas":     {"semiMajorAxisM": 1.85539e8, "eccentricity": 0.0196,  "massKg": 3.7493e19,   "retrograde": False},
    "Enceladus": {"semiMajorAxisM": 2.37948e8, "eccentricity": 0.0047,  "massKg": 1.08022e20,  "retrograde": False},
    "Tethys":    {"semiMajorAxisM": 2.94619e8, "eccentricity": 0.0001,  "massKg": 6.17449e20,  "retrograde": False},
    "Dione":     {"semiMajorAxisM": 3.77396e8, "eccentricity": 0.0022,  "massKg": 1.095452e21, "retrograde": False},
    "Rhea":      {"semiMajorAxisM": 5.27108e8, "eccentricity": 0.001258,"massKg": 2.306518e21, "retrograde": False},
    "Titan":     {"semiMajorAxisM": 1.221870e9,"eccentricity": 0.0288,  "massKg": 1.3452e23,   "retrograde": False},
    "Iapetus":   {"semiMajorAxisM": 3.560820e9,"eccentricity": 0.028612,"massKg": 1.805635e21, "retrograde": False},
}

# Fixed (not CLI-configurable) seed for each moon's starting orbital phase,
# so the moons don't all launch lined up at time zero but the layout is
# reproducible across runs and restarts.
MOON_PHASE_SEED = 20260706

MOON_COLORS = {
    "Mimas": "red", "Enceladus": "yellow", "Tethys": "cyan", "Dione": "magenta",
    "Rhea": "lime", "Titan": "gold", "Iapetus": "deepskyblue",
}


if __name__ == "__main__":
    commandLine = CommandLineArguments(animate=True,
                                        nBodies = 16000,
                                        cmass = 1.0,
                                        cycles = 1200000,
                                        rootName = "saturn-rings",
                                        restartDir = "sr/restart",
                                        vizDir = "sr/viz",
                                        restartCycle = 500,
                                        restoreCycle = None,
                                        vizCycle = 500)
    restartDir = restartDir+"-nBodies=%d"%nBodies
    vizDir = vizDir+"-nBodies=%d"%nBodies
    os.makedirs(restartDir, exist_ok=True)
    os.makedirs(vizDir, exist_ok=True)

    constants = PhysicalConstants(58232e+3,     # saturn radius in m
                                  5.683e+26,    # saturn mass in kg
                                  3.156e7       # 1 year in s
                                  )
    print("G: ",constants.G)

    cmLoc = Vector2d(0,0)

    # Saturn's dense rings (D through F) span roughly these radii in reality;
    # everything past the F ring out to the moons is essentially empty space
    # (aside from the tenuous, dust-only G/E rings).
    ringInnerRadius = 66.9e6 / constants.unitLengthMeters    # D ring inner edge
    ringOuterRadius = 140.18e6 / constants.unitLengthMeters  # F ring outer edge

    from PoissonNodeGenerator import PoissonDisk2d
    pd = PoissonDisk2d(nBodies,seed=1)

    for i in range(len(pd.positions)):
        pd.positions[i] = (pd.positions[i][0]*ringOuterRadius,pd.positions[i][1]*ringOuterRadius)
    pd.positions = [position for position in pd.positions if np.sqrt(position[0]**2 + position[1]**2) > ringInnerRadius]

    nBodies = len(pd.positions)
    myNodeList = NodeList(nBodies)

    kinematics = Kinematics2d(nodeList=myNodeList, constants=constants)
    saturnGrav = PointSourceGravity2d(nodeList=myNodeList,
                                      constants=constants,
                                      pointSourceLocation=cmLoc,
                                      pointSourceMass=cmass,
                                      pointSourceVelocity = Vector2d(0,0))

    moonGM = cmass * constants.G
    moonPhaseRng = random.Random(MOON_PHASE_SEED)

    packages = [kinematics, saturnGrav]
    moonGravs = {}
    moonMovers = []
    for name, data in SATURN_MOONS.items():
        a = data["semiMajorAxisM"] / constants.unitLengthMeters
        e = data["eccentricity"]
        mass = data["massKg"] / constants.unitMassKg
        retrograde = data["retrograde"]
        M0 = moonPhaseRng.uniform(0.0, 2.0*np.pi)   # random starting phase, not all lined up at t=0

        x0, y0 = keplerPosition(0.0, moonGM, a, e, M0, retrograde)
        grav = PointSourceGravity2d(nodeList=myNodeList,
                                    constants=constants,
                                    pointSourceLocation=Vector2d(x0, y0),
                                    pointSourceMass=mass,
                                    pointSourceVelocity=Vector2d(0,0))
        moonGravs[name] = grav
        moonMovers.append(MoonMover(grav, moonGM, a, e, M0=M0, retrograde=retrograde, workCycle=1))
        packages.append(grav)

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

    periodicWork = list(moonMovers)

    restoreIfAvailable(myNodeList, integrator, restartDir, rootName, restoreCycle)

    # Prime every moon's location at the (possibly restored) current time,
    # before the main loop runs -- otherwise the first step after a restart
    # would pull the rings toward each moon's t=0 position instead of its
    # actual current position.
    for mover in moonMovers:
        mover(integrator.cycle, integrator.time, integrator.dt)

    restartWriter = RestartWriter(myNodeList, integrator)
    def dropRestart(cycle, time, dt):
        restartWriter.write(restartFileName(restartDir, rootName, cycle))
    dropRestart.cycle = restartCycle
    periodicWork += [dropRestart]

    meshWriter = SiloDump(baseName=os.path.join(vizDir, rootName),
                            nodeList=myNodeList,
                            fieldNames=["velocity"],
                            dumpCycle=vizCycle)
    periodicWork += [meshWriter]

    print(myNodeList.numNodes)

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=100)

    if animate:
        from Animation import *
        bounds = (-12, 12, -12, 12)
        AnimateScatter(bounds, stepper=controller, positions=pos, frames=10, interval=50, size=2,
                    extra_points=lambda: [(cmLoc.x, cmLoc.y)] +
                                        [(moonGravs[name].pointSourceLocation.x, moonGravs[name].pointSourceLocation.y)
                                            for name in SATURN_MOONS],
                    extra_colors=['orange'] + [MOON_COLORS[name] for name in SATURN_MOONS], extra_size=80)
    else:
        controller.Step(cycles)