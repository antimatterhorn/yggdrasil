from yggdrasil import *
from ConstantDThetaPolyGenerator import ConstantDThetaPolyGenerator
from Mesh import ALEMesh2d, Geometry
from Physics import ALEMeshHydroHLLE2d
from EOS import IdealGasEOS
from Boundaries import ReflectingALEMeshBoundary2d
from Utilities import SiloDump
import math
import os

# Sedov point-blast on a genuinely unstructured, cylindrically-conformal
# quarter-disk ALEMesh (r>=0, z>=0, r^2+z^2<=R^2) -- a quad-dominant polar
# mesh (ConstantDThetaPolyGenerator: ntheta grows outward to hold arc-length
# spacing roughly constant, doubling at intervals via the standard
# 1-inner-segment-to-2-outer-segments triangle transition, plain quads
# elsewhere, and a triangle fan from the single center point to the
# innermost ring), not a quad grid and not a Voronoi/Delaunay triangulation.
# Via axisymmetric revolution about z, this quarter (using both the r=0
# axis and a z=0 equatorial symmetry wall) represents a true spherical
# (Sedov-Taylor d=3) point explosion -- exercising real curved-boundary
# topology together with the RZ metric/source-term/axis-boundary machinery
# all at once.

GAMMA = 1.4
SEDOV_EXPONENT_3D = 2.0 / 5.0  # R ~ t^(2/(d+2)), d=3 for a true spherical blast


def buildQuarterDiskMesh(nrings, ntheta0, R):
    gen = ConstantDThetaPolyGenerator(nrings, ntheta0, thetaMin=0.0, thetaMax=math.pi / 2)

    mesh = ALEMesh2d(Geometry.CylindricalRZ)
    for p in gen.positions:
        mesh.addNode(Vector2d(p[0] * R, p[1] * R))
    for cell in gen.cells:
        mesh.addCell(cell)
    mesh.computeFaces()
    return mesh


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
    commandLine = CommandLineArguments(nrings=32, ntheta0=5, R=3.0,
                                       siloDump = True,
                                       dumpCycle = 50,
                                       cycles=10000,
                                       rootName = "aleSedovRZ",
                                       vizDir = "viz",
                                       restartDir = "restart")

    dumpDir = rootName + "/" + "nrings=%d_ntheta0=%d_R=%3.2f" % (nrings, ntheta0, R)

    vizDir = dumpDir + "/" + vizDir
    restartDir = dumpDir + "/" + restartDir
    os.makedirs(vizDir, exist_ok=True)
    os.makedirs(restartDir, exist_ok=True)

    mesh = buildQuarterDiskMesh(nrings, ntheta0, R)
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

    # Deposit blast energy into the entire innermost ring (the ntheta0
    # triangle-fan cells around the center -- ConstantDThetaPolyGenerator
    # always builds these first, so they're cells 0..ntheta0-1). Filling
    # only a handful of cells nearest the origin by distance would cover a
    # narrow angular wedge instead of the full ring, breaking axisymmetry
    # and making the blast visibly spread angularly before advecting to the
    # next ring out.
    E_blob = 2000.0
    for i in range(ntheta0):
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
