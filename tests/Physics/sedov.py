from yggdrasil import *
import numpy as np
import matplotlib.pyplot as plt
from Physics import GridHydroHLLE2d
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
# center is compared against the analytic profile (L1 error), and total
# energy, conserved by the reflecting-walled finite-volume scheme, is checked
# as an invariant.

GAMMA = 1.4
GEOMETRY = 2                                # cylindrical symmetry of a 2D (line) blast


def totalEnergy(grid, density, velocity, energy, nx, ny):
    E = 0.0
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            v = velocity[idx]
            E += density[idx] * (energy[idx] + 0.5 * (v.x * v.x + v.y * v.y)) * grid.cellVolume(idx)
    return E


if __name__ == "__main__":
    commandLine = CommandLineArguments(cycles = 20000,
                                       nx = 100,
                                       ny = 100,
                                       dx = 1.0,
                                       dy = 1.0,
                                       tstop = 12.0,
                                       dtmin = 1e-7,
                                       plot = True)

    myGrid = Grid2d(nx, ny, dx, dy)
    myNodeList = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroHLLE2d(myNodeList, constants, eos, myGrid)
    boundary = ReflectingGridBoundary2d(grid=myGrid)
    hydro.addBoundary(boundary)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")

    x0, y0 = nx // 2, ny // 2
    sigma = 1.5
    e0 = 1000.0
    rho_ambient = 1.0
    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            r2 = (i - x0) ** 2 + (j - y0) ** 2
            energy.setValue(idx, e0 * np.exp(-r2 / (2.0 * sigma ** 2)) + 1e-3)
            density.setValue(idx, rho_ambient)

    E0 = totalEnergy(myGrid, density, velocity, energy, nx, ny)

    controller = Controller(integrator=integrator, periodicWork=[], statStep=50, tstop=tstop)
    controller.Step(cycles)
    t = controller.time

    pressure = myNodeList.getFieldDouble("pressure")
    E1 = totalEnergy(myGrid, density, velocity, energy, nx, ny)
    dE = abs(E1 - E0) / E0

    sedov = SedovSolution(GAMMA, GEOMETRY, eblast=E0, rho0=rho_ambient)
    Rshock = sedov.shockRadius(t)

    # Radial lineout from the blast center along +x -- unaffected by the
    # reflecting walls as long as the shock stays well inside the domain.
    r_sim = np.array([(i - x0) * dx for i in range(x0, nx)])
    rho_sim = np.array([density[myGrid.index(i, y0, 0)] for i in range(x0, nx)])
    p_sim = np.array([pressure[myGrid.index(i, y0, 0)] for i in range(x0, nx)])
    u_sim = np.array([velocity[myGrid.index(i, y0, 0)].x for i in range(x0, nx)])

    rho_ex, u_ex, p_ex = sedov.profile(t, r_sim)
    rho2, u2, p2 = sedov.shockJump(t)

    # L1 errors normalized by the exact post-shock jump scale (not the max of
    # the sampled arrays, which can miss the narrow near-shock spike), so the
    # large unshocked-ambient region (already near-exact) doesn't wash out the
    # error in the thin, under-resolved shocked shell that actually matters.
    err_rho = np.mean(np.abs(rho_sim - rho_ex)) / (rho2 - rho_ambient)
    err_p = np.mean(np.abs(p_sim - p_ex)) / p2
    err_u = np.mean(np.abs(u_sim - u_ex)) / u2

    print(f"Sedov at t={t:.3f} (cycle {controller.cycle}): shock R_exact={Rshock:.2f}; "
          f"L1 errors (normalized) rho={err_rho:.3e} p={err_p:.3e} u={err_u:.3e}; dE/E={dE:.2e}")
    assert Rshock < 0.45 * nx * dx, "shock reached the domain boundary; measurement invalid"
    assert err_rho < 0.15, f"density profile off: {err_rho:.3e}"
    assert err_p < 0.15, f"pressure profile off: {err_p:.3e}"
    assert err_u < 0.15, f"velocity profile off: {err_u:.3e}"
    assert dE < 2.0e-2, f"energy not conserved: dE/E={dE:.2e}"
    print("[ok] Sedov blast matches the self-similar analytic profile and conserves energy.")

    if plot:
        r_fine = np.linspace(0.0, r_sim.max(), 400)
        rho_fine, u_fine, p_fine = sedov.profile(t, r_fine)

        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        for ax, sim, fine, label in ((axes[0], rho_sim, rho_fine, "density"),
                                     (axes[1], p_sim, p_fine, "pressure"),
                                     (axes[2], u_sim, u_fine, "radial velocity")):
            ax.plot(r_sim, sim, "o", ms=3, label="simulation")
            ax.plot(r_fine, fine, "-", label="Sedov analytic")
            ax.set_xlabel("r"); ax.set_ylabel(label); ax.legend()
        fig.suptitle(f"Sedov blast (2D cylindrical) at t={t:.3f}, cycle={controller.cycle}")
        fig.tight_layout()
        plt.show()
