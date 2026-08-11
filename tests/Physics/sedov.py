from yggdrasil import *
import numpy as np
import matplotlib.pyplot as plt
from Physics import GridHydroKT2d, GridHydroHLLC2d, GridHydroHLLE2d
from Mesh import Grid2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d
from AnalyticSolutions import SedovSolution

# Sedov-Taylor point blast on a 2D Cartesian grid (cylindrically symmetric,
# i.e. a line explosion: geometry index 2), verified against the closed-form
# self-similar solution (see _sedovAnalytic.py) at a fixed physical time
# tstop -- rather than a fixed cycle count, since the analytic profile is only
# defined at a specific time, not an arbitrary cycle that shifts with
# resolution/dt. A radial lineout of density/pressure/velocity from the blast
# center is compared against the analytic profile (L1 error) for all three grid
# solvers (KT, HLLC, HLLE), and total energy, conserved by the
# reflecting-walled finite-volume scheme, is checked as an invariant.

GAMMA = 1.4
GEOMETRY = 2                                # cylindrical symmetry of a 2D (line) blast

# Profile-error tolerance per solver, normalized by the exact post-shock jump.
SOLVERS = [("KT", GridHydroKT2d, 0.15),
           ("HLLC", GridHydroHLLC2d, 0.15),
           ("HLLE", GridHydroHLLE2d, 0.15)]
DE_TOL = 2.0e-2


def totalEnergy(grid, density, velocity, energy, nx, ny):
    E = 0.0
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            v = velocity[idx]
            E += density[idx] * (energy[idx] + 0.5 * (v.x * v.x + v.y * v.y)) * grid.cellVolume(idx)
    return E


def setup(Solver, nx, ny, dx, dy, dtmin, e0, sigma, rho_ambient):
    """Build a complete Sedov problem for one solver. Every object is returned and
    kept alive by the caller: the physics package stores non-owning pointers."""
    grid = Grid2d(nx, ny, dx, dy)
    nodeList = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = Solver(nodeList, constants, eos, grid)
    boundary = ReflectingGridBoundary2d(grid=grid)
    hydro.addBoundary(boundary)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    density = nodeList.getFieldDouble("density")
    energy = nodeList.getFieldDouble("specificInternalEnergy")
    x0, y0 = nx // 2, ny // 2
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            r2 = (i - x0) ** 2 + (j - y0) ** 2
            energy.setValue(idx, e0 * np.exp(-r2 / (2.0 * sigma ** 2)) + 1e-3)
            density.setValue(idx, rho_ambient)

    return dict(grid=grid, nodeList=nodeList, constants=constants, eos=eos,
                hydro=hydro, boundary=boundary, integrator=integrator, x0=x0, y0=y0)


if __name__ == "__main__":
    commandLine = CommandLineArguments(cycles = 20000,
                                       nx = 100,
                                       ny = 100,
                                       dx = 1.0,
                                       dy = 1.0,
                                       tstop = 12.0,
                                       dtmin = 1e-7,
                                       plot = True)

    sigma = 1.5
    e0 = 1000.0
    rho_ambient = 1.0

    results = []
    failures = []
    for label, Solver, tol in SOLVERS:
        problem = setup(Solver, nx, ny, dx, dy, dtmin, e0, sigma, rho_ambient)
        grid, nodeList = problem["grid"], problem["nodeList"]
        x0, y0 = problem["x0"], problem["y0"]
        density = nodeList.getFieldDouble("density")
        energy = nodeList.getFieldDouble("specificInternalEnergy")
        velocity = nodeList.getFieldVector2d("velocity")

        E0 = totalEnergy(grid, density, velocity, energy, nx, ny)

        controller = Controller(integrator=problem["integrator"], periodicWork=[],
                                statStep=50, tstop=tstop)
        controller.Step(cycles)
        t = controller.time

        pressure = nodeList.getFieldDouble("pressure")
        E1 = totalEnergy(grid, density, velocity, energy, nx, ny)
        dE = abs(E1 - E0) / E0

        sedov = SedovSolution(GAMMA, GEOMETRY, eblast=E0, rho0=rho_ambient)
        Rshock = sedov.shockRadius(t)

        # Radial lineout from the blast center along +x -- unaffected by the
        # reflecting walls as long as the shock stays well inside the domain.
        r_sim = np.array([(i - x0) * dx for i in range(x0, nx)])
        rho_sim = np.array([density[grid.index(i, y0, 0)] for i in range(x0, nx)])
        p_sim = np.array([pressure[grid.index(i, y0, 0)] for i in range(x0, nx)])
        u_sim = np.array([velocity[grid.index(i, y0, 0)].x for i in range(x0, nx)])

        rho_ex, u_ex, p_ex = sedov.profile(t, r_sim)
        rho2, u2, p2 = sedov.shockJump(t)

        # L1 errors normalized by the exact post-shock jump scale (not the max of
        # the sampled arrays, which can miss the narrow near-shock spike), so the
        # large unshocked-ambient region (already near-exact) doesn't wash out the
        # error in the thin, under-resolved shocked shell that actually matters.
        err_rho = np.mean(np.abs(rho_sim - rho_ex)) / (rho2 - rho_ambient)
        err_p = np.mean(np.abs(p_sim - p_ex)) / p2
        err_u = np.mean(np.abs(u_sim - u_ex)) / u2

        print(f"{label} Sedov at t={t:.3f} (cycle {controller.cycle}): shock R_exact={Rshock:.2f}; "
              f"L1 errors (normalized) rho={err_rho:.3e} p={err_p:.3e} u={err_u:.3e}; "
              f"dE/E={dE:.2e}  (tol {tol:.2f})")
        if not Rshock < 0.45 * nx * dx:
            failures.append(f"{label}: shock reached the domain boundary; measurement invalid")
        for name, err in (("density", err_rho), ("pressure", err_p), ("velocity", err_u)):
            if not err < tol:
                failures.append(f"{label} {name} profile off: {err:.3e} (tol {tol:.2f})")
        if not dE < DE_TOL:
            failures.append(f"{label} energy not conserved: dE/E={dE:.2e} (tol {DE_TOL:.0e})")

        results.append((label, t, r_sim, {"density": rho_sim,
                                          "pressure": p_sim,
                                          "radial velocity": u_sim}, sedov))

    assert not failures, "\n".join(failures)
    print("[ok] Sedov blast matches the self-similar analytic profile and conserves "
          "energy for KT, HLLC, and HLLE.")

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
            ax.plot(r_fine, fine, "k-", lw=1, label="Sedov analytic")
            for label, t, r_sim, fields, _ in results:
                ax.plot(r_sim, fields[name], "o", ms=2.5, alpha=0.7, label=label)
            ax.set_xlabel("r"); ax.set_ylabel(name); ax.legend()
        fig.suptitle(f"Sedov blast (2D cylindrical) at t={t_plot:.3f}: KT vs HLLC vs HLLE")
        fig.tight_layout()
        plt.show()
