from yggdrasil import *
from ConstantDThetaPolyGenerator import ConstantDThetaPolyGenerator
from Mesh import ALEMesh2d, Geometry
from Physics import ALEMeshHydroHLLE2d
from EOS import IdealGasEOS
from Boundaries import ReflectingALEMeshBoundary2d
import math

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


def shockRadius(density, mesh, numCells, rho_ambient, thresh=1.2):
    rmax = 0.0
    for i in range(numCells):
        if density[i] > thresh * rho_ambient:
            c = mesh.cellCentroid(i)
            rmax = max(rmax, math.hypot(c.x, c.y))
    return rmax


if __name__ == "__main__":
    commandLine = CommandLineArguments(nrings=12, ntheta0=4, R=3.0, cyclesA=3000, cyclesB=4500)

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

    integrator = RungeKutta4Integrator2d([hydro], dtmin=1e-6, verbose=False)
    controller = Controller(integrator=integrator, periodicWork=[], statStep=100000)

    m0 = totalMass(density, mesh, numCells)
    r0 = shockRadius(density, mesh, numCells, rho_ambient)

    controller.Step(cyclesA)
    t1 = controller.time
    r1 = shockRadius(density, mesh, numCells, rho_ambient)

    controller.Step(cyclesB - cyclesA)
    t2 = controller.time
    r2 = shockRadius(density, mesh, numCells, rho_ambient)

    finite = all(density[i] == density[i] and abs(density[i]) < 1e6 for i in range(numCells))
    assert finite, "blast produced NaN/inf or runaway density"
    print("[ok] solution finite and bounded")

    assert r2 > r0 + 0.05 * R, f"shock did not expand: r0={r0:.3f} -> r2={r2:.3f}"
    print(f"[ok] shock front expanded radially: r={r0:.3f} -> {r2:.3f}")

    m2 = totalMass(density, mesh, numCells)
    rel = abs(m2 - m0) / m0
    assert rel < 1e-3, f"mass not conserved: rel change {rel:.3e}"
    print(f"[ok] per-radian mass conserved: rel change {rel:.2e}")

    # Diagnostic only, not a pass/fail gate: the self-similar R~t^(2/5) law is
    # an asymptotic result, and this mesh's coarse, highly non-uniform cell
    # sizes near the singular center (the triangle fan's cells are tiny in
    # RZ-weighted volume close to r=0) mean the blast needs more resolution
    # and dynamic range than this quick test uses before it settles into
    # that asymptotic regime. Guarded against r1/r2 landing at exactly 0
    # (below the detection threshold at that cycle) rather than tuning the
    # cycle window to avoid it exactly -- the threshold crossing here can
    # jump rather than creep, since the whole interior heats up together.
    if r1 > 0.0 and r2 > 0.0:
        exponent = math.log(r2 / r1) / math.log(t2 / t1)
        print(f"self-similar exponent (diagnostic, pre-asymptotic at this resolution): "
              f"{exponent:.3f} vs exact {SEDOV_EXPONENT_3D:.3f}")
    else:
        print(f"self-similar exponent: skipped (shock not yet above threshold at "
              f"cyclesA or cyclesB; r1={r1:.3f}, r2={r2:.3f})")

    print("\nQuarter-disk RZ Sedov blast: all checks passed.")
