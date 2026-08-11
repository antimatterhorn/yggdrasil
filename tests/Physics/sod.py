from yggdrasil import *
import numpy as np
import matplotlib.pyplot as plt
from Animation import *
from Physics import GridHydroKT2d, GridHydroHLLC2d, GridHydroHLLE2d
from Mesh import Grid2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d
from AnalyticSolutions import SodSolution

# Sod shock tube, verified against the exact Riemann solution (AnalyticSolutions.py),
# for all three grid solvers (KT, HLLC, HLLE).
# The classic initial state (rho,u,p) = (1,0,1) | (0.125,0,0.1) with gamma=1.4
# develops a left rarefaction, a contact, and a right shock; the exact
# self-similar solution is sampled at the achieved time and compared to each
# numerical profile (L1 error), and all three are plotted together when plot=True.

GAMMA = 1.4
RHOL, PL, RHOR, PR = 1.0, 1.0, 0.125, 0.1

# HLLC/HLLE carry a slightly larger L1 error than KT on this problem, mostly at
# the contact, so each solver gets its own tolerance.
SOLVERS = [("KT", GridHydroKT2d, 2.0e-2),
           ("HLLC", GridHydroHLLC2d, 3.0e-2),
           ("HLLE", GridHydroHLLE2d, 3.0e-2)]


def setup(Solver, nx, ny, dx, dy, dtmin):
    """Build a complete Sod problem for one solver. Every object is returned and
    kept alive by the caller: the physics package stores non-owning pointers."""
    grid = Grid2d(nx, ny, dx, dy)
    nodeList = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = Solver(nodeList, constants, eos, grid)
    boundary = ReflectingGridBoundary2d(grid=grid)
    hydro.addBoundary(boundary)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    x0 = 0.5 * nx * dx
    density = nodeList.getFieldDouble("density")
    energy = nodeList.getFieldDouble("specificInternalEnergy")
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            if (i + 0.5) * dx < x0:
                density.setValue(idx, RHOL); energy.setValue(idx, PL / ((GAMMA - 1) * RHOL))
            else:
                density.setValue(idx, RHOR); energy.setValue(idx, PR / ((GAMMA - 1) * RHOR))

    return dict(grid=grid, nodeList=nodeList, constants=constants, eos=eos,
                hydro=hydro, boundary=boundary, integrator=integrator, x0=x0)


def lineout(problem, nx, ny, dx):
    """Mid-plane (x, rho, p, ux) profile of a finished run."""
    grid, nodeList = problem["grid"], problem["nodeList"]
    density = nodeList.getFieldDouble("density")
    pressure = nodeList.getFieldDouble("pressure")
    velocity = nodeList.getFieldVector2d("velocity")
    jmid = ny // 2
    x = np.array([(i + 0.5) * dx for i in range(nx)])
    rho = np.array([density[grid.index(i, jmid, 0)] for i in range(nx)])
    p = np.array([pressure[grid.index(i, jmid, 0)] for i in range(nx)])
    u = np.array([velocity[grid.index(i, jmid, 0)].x for i in range(nx)])
    return x, rho, p, u


if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       animateSolver = "KT",
                                       cycles = 20000,
                                       nx = 200,
                                       ny = 10,
                                       dx = 0.005,
                                       dy = 0.005,
                                       tstop = 0.2,
                                       dtmin = 1e-7,
                                       plot = True)

    if animate:
        byName = dict((label, cls) for label, cls, _ in SOLVERS)
        assert animateSolver in byName, f"--animateSolver must be one of {list(byName)}"
        Solver = byName[animateSolver]
        problem = setup(Solver, nx, ny, dx, dy, dtmin)
        hydro = problem["hydro"]
        controller = Controller(integrator=problem["integrator"], periodicWork=[],
                                statStep=200, tstop=tstop)
        title = MakeTitle(controller, "time", "time")
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d, stepper=controller.Step,
                                                title=title, fieldName="density")
        AnimateGrid2d((nx, ny), update_method, extremis=[0, 1.1], frames=cycles, cmap="plasma")
    else:
        sod = SodSolution(GAMMA, RHOL, 0.0, PL, RHOR, 0.0, PR)
        results = []
        failures = []
        for label, Solver, tol in SOLVERS:
            problem = setup(Solver, nx, ny, dx, dy, dtmin)
            controller = Controller(integrator=problem["integrator"], periodicWork=[],
                                    statStep=200, tstop=tstop)
            controller.Step(cycles)
            t = controller.time
            x_sim, rho_sim, p_sim, u_sim = lineout(problem, nx, ny, dx)

            rho_ex, u_ex, p_ex = sod.profile(t, x_sim, x0=problem["x0"])
            err_rho = np.mean(np.abs(rho_sim - rho_ex))
            err_u = np.mean(np.abs(u_sim - u_ex))
            err_p = np.mean(np.abs(p_sim - p_ex))

            print(f"{label} Sod at t={t:.4f} (cycle {controller.cycle}): "
                  f"L1 errors  rho={err_rho:.4e}  u={err_u:.4e}  p={err_p:.4e}  (tol {tol:.0e})")
            for name, err in (("density", err_rho), ("velocity", err_u), ("pressure", err_p)):
                if not err < tol:
                    failures.append(f"{label} {name} L1 error too large: {err:.4e} (tol {tol:.0e})")
            results.append((label, t, x_sim, {"density": rho_sim,
                                              "pressure": p_sim,
                                              "velocity": u_sim}))

        assert not failures, "\n".join(failures)
        print("[ok] Sod shock tube matches the exact Riemann solution for KT, HLLC, and HLLE.")

        if plot:
            # Each solver overshoots tstop by at most one dt, so their achieved
            # times differ negligibly; the exact curve is drawn at the first one's.
            t_plot = results[0][1]
            x_fine = np.linspace(results[0][2].min(), results[0][2].max(), 400)
            rho_fine, u_fine, p_fine = sod.profile(t_plot, x_fine, x0=0.5 * nx * dx)

            fig, axes = plt.subplots(1, 3, figsize=(14, 4))
            for ax, fine, name in ((axes[0], rho_fine, "density"),
                                   (axes[1], p_fine, "pressure"),
                                   (axes[2], u_fine, "velocity")):
                ax.plot(x_fine, fine, "k-", lw=1, label="exact Riemann")
                for label, t, x_sim, fields in results:
                    ax.plot(x_sim, fields[name], "o", ms=2.5, alpha=0.7, label=label)
                ax.set_xlabel("x"); ax.set_ylabel(name); ax.legend()
            fig.suptitle(f"Sod shock tube at t={t_plot:.4f}: KT vs HLLC vs HLLE")
            fig.tight_layout()
            plt.show()
