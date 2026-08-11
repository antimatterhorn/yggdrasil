from yggdrasil import *
import math
import numpy as np
import matplotlib.pyplot as plt
from Physics import GridHydroKT2d, GridHydroHLLC2d, GridHydroHLLE2d
from Mesh import Grid2d, Geometry
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d
from AnalyticSolutions import SedovSolution

# Sedov-like point blast at the (r=0, z=0) corner of an axisymmetric (r,z)
# grid, run with all three grid solvers (KT, HLLC, HLLE). Together, the
# auto-installed r=0 axis boundary and a reflecting z=0 equatorial wall
# represent one hemisphere of a true spherical (Sedov-Taylor geometry index 3)
# point explosion -- revolving the domain about z sweeps the full azimuthal
# angle, so only the z>=0 restriction is an actual symmetry cut. This exercises
# the full RZ path (radius-weighted volumes/areas, the p/r source, and the axis
# boundary) on a genuinely radial problem, verified against the closed-form
# self-similar solution (see _sedovAnalytic.py) at a fixed physical time tstop,
# sampled along the z=0 row where cell radius equals the true spherical radius.

GAMMA = 1.4
GEOMETRY = 3                                # spherical: axis + equatorial wall together

# Profile-error tolerance per solver, normalized by the exact post-shock jump.
SOLVERS = [("KT", GridHydroKT2d, 0.15),
           ("HLLC", GridHydroHLLC2d, 0.15),
           ("HLLE", GridHydroHLLE2d, 0.15)]


def interior_ids(grid, nr, nz):
    return [grid.index(i, j, 0) for j in range(nz) for i in range(nr)
            if not grid.onBoundary(grid.index(i, j, 0))]


def total_mass(density, grid, ids):
    return sum(density[idx] * grid.cellVolume(idx) for idx in ids)


def total_energy_halfsphere(density, energy, velocity, grid, ids):
    """Energy summed over the simulated (z>=0, per-radian) domain. Multiply by
    4*pi (2*pi azimuthal revolution x2 for the z<0 mirror half) to get the
    true total blast energy of the full sphere the RZ domain represents."""
    E = 0.0
    for idx in ids:
        v = velocity[idx]
        E += density[idx] * (energy[idx] + 0.5 * (v.x * v.x + v.y * v.y)) * grid.cellVolume(idx)
    return E


def shock_radius(density, grid, nr, nz, rho_ambient):
    rmax = 0.0
    for j in range(nz):
        for i in range(nr):
            idx = grid.index(i, j, 0)
            if density[idx] > 1.5 * rho_ambient:
                rmax = max(rmax, grid.cellRadius(idx))
    return rmax


def setup(Solver, nr, nz, dr, dtmin, rho_ambient, blob, E_blob):
    """Build a complete RZ Sedov problem for one solver. Every object is returned
    and kept alive by the caller: the physics package stores non-owning pointers."""
    grid = Grid2d(nr, nz, dr, dr, Geometry.CylindricalRZ)
    nodeList = NodeList(nr * nz)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = Solver(nodeList, constants, eos, grid)
    box = ReflectingGridBoundary2d(grid=grid)
    box.setFaces(["right", "top", "bottom"])   # 'left' (r=0 axis) auto-installed
    hydro.addBoundary(box)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    density = nodeList.getFieldDouble("density")
    energy = nodeList.getFieldDouble("specificInternalEnergy")
    for i in range(nr * nz):
        density.setValue(i, rho_ambient)
        energy.setValue(i, 1e-3)               # cold ambient

    # Deposit energy in a blob at the axis/base corner.
    for j in range(blob):
        for i in range(blob):
            energy.setValue(grid.index(i, j, 0), E_blob)

    return dict(grid=grid, nodeList=nodeList, constants=constants, eos=eos,
                hydro=hydro, boundary=box, integrator=integrator)


if __name__ == "__main__":
    commandLine = CommandLineArguments(nr = 60,
                                       nz = 60,
                                       dr = 0.02,
                                       cycles = 20000,
                                       tstop = 0.3,
                                       dtmin = 1e-6,
                                       plot = True)

    rho_ambient = 1.0
    blob = 2
    E_blob = 400.0

    results = []
    failures = []
    for label, Solver, tol in SOLVERS:
        problem = setup(Solver, nr, nz, dr, dtmin, rho_ambient, blob, E_blob)
        grid, nodes = problem["grid"], problem["nodeList"]
        density = nodes.getFieldDouble("density")
        energy = nodes.getFieldDouble("specificInternalEnergy")
        velocity = nodes.getFieldVector2d("velocity")

        ids = interior_ids(grid, nr, nz)
        m0 = total_mass(density, grid, ids)
        r0 = shock_radius(density, grid, nr, nz, rho_ambient)
        Eblast = 4.0 * math.pi * total_energy_halfsphere(density, energy, velocity, grid, ids)

        controller = Controller(integrator=problem["integrator"], periodicWork=[],
                                statStep=50, tstop=tstop)
        controller.Step(cycles)
        t = controller.time

        # (i) finite and bounded
        if not all(density[i] == density[i] and abs(density[i]) < 1e6
                   for i in range(nr * nz)):
            failures.append(f"{label}: blast produced NaN/inf or runaway density")

        # (ii) shock expanded outward
        r1 = shock_radius(density, grid, nr, nz, rho_ambient)
        if not r1 > r0 + 5 * dr:
            failures.append(f"{label}: shock did not expand: r0={r0:.3f} -> r1={r1:.3f}")

        # (iii) per-radian mass conserved
        m1 = total_mass(density, grid, ids)
        rel = abs(m1 - m0) / m0
        if not rel < 1e-3:
            failures.append(f"{label}: mass not conserved: rel change {rel:.3e}")

        # (iv) matches the self-similar spherical Sedov profile, sampled along the
        # z=0 row where cell radius equals the true spherical radius.
        sedov = SedovSolution(GAMMA, GEOMETRY, eblast=Eblast, rho0=rho_ambient)
        Rshock = sedov.shockRadius(t)

        pressure = nodes.getFieldDouble("pressure")
        r_sim = np.array([grid.cellRadius(grid.index(i, 0, 0)) for i in range(nr)])
        rho_sim = np.array([density[grid.index(i, 0, 0)] for i in range(nr)])
        p_sim = np.array([pressure[grid.index(i, 0, 0)] for i in range(nr)])
        u_sim = np.array([velocity[grid.index(i, 0, 0)].x for i in range(nr)])

        rho_ex, u_ex, p_ex = sedov.profile(t, r_sim)
        rho2, u2, p2 = sedov.shockJump(t)
        err_rho = np.mean(np.abs(rho_sim - rho_ex)) / (rho2 - rho_ambient)
        err_p = np.mean(np.abs(p_sim - p_ex)) / p2
        err_u = np.mean(np.abs(u_sim - u_ex)) / u2

        print(f"{label} Sedov RZ at t={t:.4f} (cycle {controller.cycle}): "
              f"shock R_exact={Rshock:.3f} (measured {r1:.3f}); mass drift {rel:.2e}; "
              f"L1 errors (normalized) rho={err_rho:.3e} p={err_p:.3e} u={err_u:.3e} "
              f"(tol {tol:.2f})")
        if not Rshock < 0.9 * nr * dr:
            failures.append(f"{label}: shock reached the domain boundary; measurement invalid")
        for name, err in (("density", err_rho), ("pressure", err_p), ("velocity", err_u)):
            if not err < tol:
                failures.append(f"{label} {name} profile off: {err:.3e} (tol {tol:.2f})")

        results.append((label, t, r_sim, {"density": rho_sim,
                                          "pressure": p_sim,
                                          "radial velocity": u_sim}, sedov))

    assert not failures, "\n".join(failures)
    print("\nCylindrical Sedov blast: all checks passed for KT, HLLC, and HLLE.")

    if plot:
        # Each solver overshoots tstop by at most one dt, so their achieved times
        # and blast energies differ negligibly; the analytic curve uses the first's.
        label0, t_plot, r0, _, sedov0 = results[0]
        r_fine = np.linspace(0.0, r0.max(), 400)
        rho_fine, u_fine, p_fine = sedov0.profile(t_plot, r_fine)

        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        for ax, fine, name in ((axes[0], rho_fine, "density"),
                               (axes[1], p_fine, "pressure"),
                               (axes[2], u_fine, "radial velocity")):
            ax.plot(r_fine, fine, "k-", lw=1, label="Sedov analytic (d=3)")
            for label, t, r_sim, fields, _ in results:
                ax.plot(r_sim, fields[name], "o", ms=2.5, alpha=0.7, label=f"{label} (z=0 row)")
            ax.set_xlabel("r"); ax.set_ylabel(name); ax.legend()
        fig.suptitle(f"Cylindrical (r,z) Sedov blast at t={t_plot:.4f}: KT vs HLLC vs HLLE")
        fig.tight_layout()
        plt.show()
