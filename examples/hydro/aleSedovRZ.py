from yggdrasil import *
from ConstantDThetaPolyGenerator import ConstantDThetaPolyDisk2d
from Mesh import ALEMesh2d, Geometry
from Physics import ALEMeshHydroHLLE2d
from EOS import IdealGasEOS
from Boundaries import ReflectingALEMeshBoundary2d
from Utilities import SiloDump
import math
import os

# Sedov point-blast on a genuinely unstructured, cylindrically-conformal
# quarter-disk ALEMesh (r>=0, z>=0, r^2+z^2<=R^2) -- a quad-dominant polar
# mesh (ConstantDThetaPolyDisk2d: Delaunay-triangulate ConstantDThetaDisk2d's
# seed points, then greedily merge adjacent triangle pairs into quads
# wherever the result is convex and close enough to square, ranked
# best-quality-first; whatever can't cleanly merge is left as a triangle),
# not a quad grid and not a fixed doubling-transition scheme. Via
# axisymmetric revolution about z, this quarter (using both the r=0 axis
# and a z=0 equatorial symmetry wall) represents a true spherical
# (Sedov-Taylor d=3) point explosion -- exercising real curved-boundary
# topology together with the RZ metric/source-term/axis-boundary machinery
# all at once.

GAMMA = 1.4
SEDOV_EXPONENT_3D = 2.0 / 5.0  # R ~ t^(2/(d+2)), d=3 for a true spherical blast


def buildQuarterDiskMesh(npoints, R, angleTol=25.0):
    gen = ConstantDThetaPolyDisk2d(npoints, thetaMin=0.0, thetaMax=math.pi / 2, angleTol=angleTol)

    mesh = ALEMesh2d(Geometry.CylindricalRZ)
    for p in gen.positions:
        mesh.addNode(Vector2d(p[0] * R, p[1] * R))
    for cell in gen.cells:
        mesh.addCell(cell)
    mesh.computeFaces()

    # Cells touching the origin node (always node 0 -- ConstantDThetaDisk2d
    # adds the center first) -- however many the triangulation happens to
    # put there, found by topology rather than assumed from cell ordering
    # or a distance cutoff, since this mesh has no ring structure to exploit.
    originCells = [i for i, cell in enumerate(gen.cells) if 0 in cell]
    return mesh, originCells


def wallFaces(mesh, R):
    # Every boundary face except the r=0 axis, which ALEMeshHydroBase
    # auto-installs a reflecting condition on for CylindricalRZ geometry.
    tiny = 1e-6 * R
    faces = mesh.getFaces()
    return [idx for idx, f in enumerate(faces) if f.isBoundary and f.centroid.x >= tiny]


def totalMass(density, mesh, numCells):
    return sum(density[i] * mesh.cellVolume(i) for i in range(numCells))


def shockRadius(density, mesh, numCells, rho_ambient):
    rmax = 0.0
    for i in range(numCells):
        if density[i] > 1.5 * rho_ambient:
            c = mesh.cellCentroid(i)
            rmax = max(rmax, math.hypot(c.x, c.y))
    return rmax


if __name__ == "__main__":
    commandLine = CommandLineArguments(npoints=2000, R=3.0,
                                       siloDump = True,
                                       dumpCycle = 50,
                                       cycles=10000,
                                       rootName = "aleSedovRZ",
                                       vizDir = "viz",
                                       restartDir = "restart")

    dumpDir = rootName + "/" + "npoints=%d_R=%3.2f" % (npoints, R)

    vizDir = dumpDir + "/" + vizDir
    restartDir = dumpDir + "/" + restartDir
    os.makedirs(vizDir, exist_ok=True)
    os.makedirs(restartDir, exist_ok=True)

    mesh, originCells = buildQuarterDiskMesh(npoints, R)
    numCells = len(mesh.getConnectivityMap())
    print("quarter-disk ALEMesh cells:", numCells)

    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    nodes = NodeList(numCells)
    hydro = ALEMeshHydroHLLE2d(nodes, constants, eos, mesh)

    box = ReflectingALEMeshBoundary2d(mesh)
    box.setFaces(wallFaces(mesh, R))
    hydro.addBoundary(box)

    rho_ambient = 1.0
    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    for i in range(numCells):
        density.setValue(i, rho_ambient)
        energy.setValue(i, 1e-3)

    # Deposit blast energy into every cell touching the origin (all of
    # originCells, found by topology in buildQuarterDiskMesh). Filling only
    # a handful of cells nearest the origin by distance -- or assuming a
    # fixed cell-index range -- would risk covering a narrow angular wedge
    # instead of the whole region around the origin, breaking axisymmetry.
    E_blob = 2000.0
    for i in originCells:
        energy.setValue(i, E_blob)

    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName=os.path.join(vizDir, rootName),
                              nodeList=nodes,
                              fieldNames=["density", "specificInternalEnergy",
                                          "pressure", "velocity"],
                              dumpCycle=dumpCycle,
                              mesh=mesh)
        periodicWork += [meshWriter]

    integrator = RungeKutta4Integrator2d([hydro], dtmin=1e-6, verbose=False)
    controller = Controller(integrator=integrator, periodicWork=periodicWork, statStep=100)

    m0 = totalMass(density, mesh, numCells)
    r0 = shockRadius(density, mesh, numCells, rho_ambient)

    controller.Step(cycles)
