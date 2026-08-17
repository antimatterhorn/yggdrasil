from yggdrasil import *
import random
from Physics import CellFlowPhysics2d, Kinematics2d
from Boundaries import BoxCollider2d
from RandomNodeGenerator import RandomNodeGenerator2d
from Animation import AnimateScatter

# CellFlow / particle-life demo: a random, generally-asymmetric per-type
# force matrix, animated live with Animation.AnimateScatter instead of
# writing Silo output.

if __name__ == "__main__":
    commandLine = CommandLineArguments(numNodes=400,
                                       numTypes=5,
                                       baseRadius=0.12,
                                       repulsion=1.0,
                                       attraction=1.0,
                                       k=6.0,
                                       forceMultiplier=10.0,
                                       damping=1.5,
                                       dtmin=0.01,
                                       seed=1,
                                       statStep=1000)

    random.seed(seed)

    constants = MKS()

    bounds = [[-1, -1], [1, 1]]
    dist = RandomNodeGenerator2d(numNodes, bounds=bounds).positions
    nodes = NodeList(numNodes)

    # random, generally-asymmetric interaction matrix -- this asymmetry is
    # what produces CellFlow's chasing/orbiting clusters rather than simple
    # clumping (a symmetric matrix would reduce to attraction-only physics).
    forceTable = [random.uniform(-1.0, 1.0) for _ in range(numTypes*numTypes)]
    radiusByType = [random.uniform(0.6, 1.4) for _ in range(numTypes)]

    cellflow = CellFlowPhysics2d(nodes, constants,
                                 numTypes=numTypes,
                                 forceTable=forceTable,
                                 radiusByType=radiusByType,
                                 baseRadius=baseRadius,
                                 repulsion=repulsion,
                                 attraction=attraction,
                                 k=k,
                                 forceMultiplier=forceMultiplier,
                                 damping=damping)
    # Kinematics (not Kinetics): Kinetics' own elastic-collision CFL reads the
    # "radius" field too, and would re-clamp dt back down to radius/velocity
    # every step
    kinematics = Kinematics2d(nodes, constants)
    nodes.insertFieldDouble("radius")

    margin = 0.1
    xlo, ylo = bounds[0]
    xhi, yhi = bounds[1]
    walls = [
        BoxCollider2d(position1=Vector2d(xlo - margin, ylo - margin), position2=Vector2d(xlo, yhi + margin), elasticity=0.9),
        BoxCollider2d(position1=Vector2d(xhi, ylo - margin), position2=Vector2d(xhi + margin, yhi + margin), elasticity=0.9),
        BoxCollider2d(position1=Vector2d(xlo - margin, ylo - margin), position2=Vector2d(xhi + margin, ylo), elasticity=0.9),
        BoxCollider2d(position1=Vector2d(xlo - margin, yhi), position2=Vector2d(xhi + margin, yhi + margin), elasticity=0.9),
    ]
    for wall in walls:
        kinematics.addBoundary(wall)

    packages = [kinematics, cellflow]

    position = nodes.getFieldVector2d("position")
    cellType = nodes.getFieldInt("cellType")
    radius = nodes.getFieldDouble("radius")
    for i in range(numNodes):
        position.setValue(i, Vector2d(dist[i][0], dist[i][1]))
        cellType.setValue(i, random.randrange(numTypes))
        radius.setValue(i, 0.01)

    # CellFlowPhysics's own EstimateTimestep only binds once particles are
    # already moving fast enough to threaten skipping past a neighbor;
    # Kinematics contributes no opinion, so dtmin is the effective, fixed
    # per-step time unit -- raise it to get visible motion per animation frame.
    integrator = RungeKutta2Integrator2d(packages=packages, dtmin=dtmin, verbose=False)
    controller = Controller(integrator=integrator, periodicWork=[], statStep=statStep)

    AnimateScatter((-1.1, 1.1, -1.1, 1.1), stepper=controller, positions=position,
                   get_color_field=lambda i: cellType[i],
                   cmap="tab10",
                   color_limits=(0, numTypes-1),
                   size=15,
                   background="black")
