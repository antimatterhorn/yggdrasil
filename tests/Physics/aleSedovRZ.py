from yggdrasil import *
from ConstantDThetaPolyGenerator import ConstantDThetaPolyDisk2d
from Mesh import ALEMesh2d, Geometry
from Physics import ALEMeshHydroHLLE2d
from EOS import IdealGasEOS
from Boundaries import ReflectingALEMeshBoundary2d
import math
import numpy as np
import matplotlib.pyplot as plt
from AnalyticSolutions import SedovSolution

# Sedov point-blast on a genuinely unstructured, cylindrically-conformal
# quarter-disk ALEMesh (r>=0, z>=0, r^2+z^2<=R^2) -- a quad-dominant polar
# mesh (ConstantDThetaPolyDisk2d: Delaunay-triangulate ConstantDThetaDisk2d's
# seed points, then greedily merge adjacent triangle pairs into quads
# wherever the result is convex and close enough to square, ranked
# best-quality-first; whatever can't cleanly merge is left as a triangle),
# not a quad grid and not a fixed doubling-transition scheme. Via
# axisymmetric revolution about z, this quarter (using both the r=0 axis
# and a z=0 equatorial symmetry wall) represents a true spherical
# (Sedov-Taylor geometry index 3) point explosion -- exercising real
# curved-boundary topology together with the RZ metric/source-term/axis-
# boundary machinery all at once. Run to a fixed physical time tstop (not a
# cycle count, per _sedovAnalytic.py) and checked/plotted against the same
# closed-form self-similar solution as sedovRZ.py, radially binned since this
# mesh has no fixed row/lineout to sample the way a Grid does.

GAMMA = 1.4
GEOMETRY = 3                                 # spherical: axis + equatorial wall together


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


def totalEnergyHalfSphere(density, energy, velocity, mesh, numCells):
    """Energy summed over the simulated (z>=0, per-radian) quarter-disk. Multiply
    by 4*pi (2*pi azimuthal revolution x2 for the z<0 mirror half) to get the
    true total blast energy of the full sphere this domain represents -- same
    bookkeeping as sedovRZ.py's Grid version, since ALEMesh's cellVolume for
    CylindricalRZ is "per-radian" too (src/Mesh/CLAUDE.md)."""
    E = 0.0
    for i in range(numCells):
        v = velocity[i]
        E += density[i] * (energy[i] + 0.5 * (v.x * v.x + v.y * v.y)) * mesh.cellVolume(i)
    return E


def shockRadius(density, mesh, numCells, rho_ambient, thresh=1.2):
    rmax = 0.0
    for i in range(numCells):
        if density[i] > thresh * rho_ambient:
            c = mesh.cellCentroid(i)
            rmax = max(rmax, math.hypot(c.x, c.y))
    return rmax


def radialProfile(values, mesh, numCells, rmax, nbins=25):
    """Volume-weighted average of a per-cell scalar over radius bins spanning
    [0, rmax] -- the unstructured-mesh equivalent of a Grid lineout. Returns
    (bin-center radii, values) for bins containing at least one cell."""
    edges = np.linspace(0.0, rmax, nbins + 1)
    sums = np.zeros(nbins)
    weights = np.zeros(nbins)
    for i in range(numCells):
        c = mesh.cellCentroid(i)
        r = math.hypot(c.x, c.y)
        if r >= rmax:
            continue
        b = min(int(r / rmax * nbins), nbins - 1)
        w = mesh.cellVolume(i)
        sums[b] += values[i] * w
        weights[b] += w
    centers = 0.5 * (edges[:-1] + edges[1:])
    valid = weights > 0
    return centers[valid], sums[valid] / weights[valid]


def radialVelocity(velocity, mesh, numCells):
    vr = np.zeros(numCells)
    for i in range(numCells):
        c = mesh.cellCentroid(i)
        r = math.hypot(c.x, c.y)
        if r > 1e-12:
            vr[i] = (velocity[i].x * c.x + velocity[i].y * c.y) / r
    return vr


if __name__ == "__main__":
    commandLine = CommandLineArguments(npoints = 400,
                                       R = 3.0,
                                       cycles = 20000,
                                       tstop = 0.2,
                                       dtmin = 1e-6,
                                       plot = True)

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

    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    velocity = nodes.getFieldVector2d("velocity")
    m0 = totalMass(density, mesh, numCells)
    r0 = shockRadius(density, mesh, numCells, rho_ambient)
    Eblast = 4.0 * math.pi * totalEnergyHalfSphere(density, energy, velocity, mesh, numCells)

    controller = Controller(integrator=integrator, periodicWork=[], statStep=50, tstop=tstop)
    controller.Step(cycles)
    t = controller.time

    finite = all(density[i] == density[i] and abs(density[i]) < 1e6 for i in range(numCells))
    assert finite, "blast produced NaN/inf or runaway density"
    print("[ok] solution finite and bounded")

    r1 = shockRadius(density, mesh, numCells, rho_ambient)
    assert r1 > r0 + 0.05 * R, f"shock did not expand: r0={r0:.3f} -> r1={r1:.3f}"
    print(f"[ok] shock front expanded radially: r={r0:.3f} -> {r1:.3f}")

    m1 = totalMass(density, mesh, numCells)
    rel = abs(m1 - m0) / m0
    assert rel < 1e-3, f"mass not conserved: rel change {rel:.3e}"
    print(f"[ok] per-radian mass conserved: rel change {rel:.2e}")

    # Matches the self-similar spherical Sedov profile, radially binned across
    # the whole unstructured mesh (volume-weighted per bin).
    sedov = SedovSolution(GAMMA, GEOMETRY, eblast=Eblast, rho0=rho_ambient)
    Rshock = sedov.shockRadius(t)

    pressure = nodes.getFieldDouble("pressure")
    vr = radialVelocity(velocity, mesh, numCells)
    rmax = 0.9 * R
    r_bins, rho_sim = radialProfile(np.array([density[i] for i in range(numCells)]), mesh, numCells, rmax)
    _, p_sim = radialProfile(np.array([pressure[i] for i in range(numCells)]), mesh, numCells, rmax)
    _, u_sim = radialProfile(vr, mesh, numCells, rmax)

    rho_ex, u_ex, p_ex = sedov.profile(t, r_bins)
    rho2, u2, p2 = sedov.shockJump(t)
    err_rho = np.mean(np.abs(rho_sim - rho_ex)) / (rho2 - rho_ambient)
    err_p = np.mean(np.abs(p_sim - p_ex)) / p2
    err_u = np.mean(np.abs(u_sim - u_ex)) / u2

    print(f"Sedov ALE RZ at t={t:.4f} (cycle {controller.cycle}): shock R_exact={Rshock:.3f} "
          f"(measured {r1:.3f}); L1 errors (normalized) rho={err_rho:.3e} p={err_p:.3e} "
          f"u={err_u:.3e}")
    assert Rshock < 0.9 * R, "shock reached the domain boundary; measurement invalid"
    assert err_rho < 0.2, f"density profile off: {err_rho:.3e}"
    assert err_p < 0.2, f"pressure profile off: {err_p:.3e}"
    assert err_u < 0.2, f"velocity profile off: {err_u:.3e}"
    print("[ok] matches the self-similar spherical Sedov profile")

    print("\nQuarter-disk RZ Sedov blast: all checks passed.")

    if plot:
        r_fine = np.linspace(0.0, rmax, 400)
        rho_fine, u_fine, p_fine = sedov.profile(t, r_fine)

        density_all = np.array([density[i] for i in range(numCells)])
        pressure_all = np.array([pressure[i] for i in range(numCells)])
        r_all = np.array([math.hypot(*(lambda c: (c.x, c.y))(mesh.cellCentroid(i))) for i in range(numCells)])

        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        for ax, cellvals, binvals, fine, label in (
                (axes[0], density_all, rho_sim, rho_fine, "density"),
                (axes[1], pressure_all, p_sim, p_fine, "pressure"),
                (axes[2], vr, u_sim, u_fine, "radial velocity")):
            ax.plot(r_all, cellvals, ".", ms=3, color="lightgray", label="cells")
            ax.plot(r_bins, binvals, "o", ms=5, label="radial bin average")
            ax.plot(r_fine, fine, "-", label="Sedov analytic (d=3)")
            ax.set_xlabel("r"); ax.set_ylabel(label); ax.legend()
        fig.suptitle(f"Quarter-disk ALE RZ Sedov blast at t={t:.4f}, cycle={controller.cycle}")
        fig.tight_layout()
        plt.show()
