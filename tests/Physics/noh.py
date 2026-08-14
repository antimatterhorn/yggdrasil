from yggdrasil import *
import numpy as np
import matplotlib.pyplot as plt
from Animation import *
from Physics import GridHydroKT2d, GridHydroHLLC2d, GridHydroHLLE2d
from Mesh import Grid2d
from EOS import IdealGasEOS
from Boundaries import OutflowGridBoundary2d
from AnalyticSolutions import NohSolution

# Planar Noh problem, verified against its exact solution (AnalyticSolutions.py)
# for all three grid solvers (KT, HLLC, HLLE).
# Two cold streams (rho=1, p~0) collide head-on at x0 with speed v0, driving a
# pair of shocks outward from the stagnation point. For gamma=5/3, v0=1 the
# exact solution is:
#   * shocked region |x-x0| < D t : rho = rho0*(g+1)/(g-1) = 4, v = 0,
#                                    p = (g-1)*rho2*v0^2/2 = 4/3
#   * unshocked region             : rho = rho0, v = +/- v0 (toward x0), p ~ 0
# with shock speed D = v0*(g-1)/2 = 1/3. The stagnation point suffers the
# well-known Noh "wall heating" over a few cells, so the plateau is checked away
# from x0.

GAMMA = 5.0 / 3.0
V0 = 1.0
RHO0 = 1.0
P0 = 1.0e-6
GEOMETRY = 1                                # planar: two independent 1D streams

SOLVERS = [("KT", GridHydroKT2d),
           ("HLLC", GridHydroHLLC2d),
           ("HLLE", GridHydroHLLE2d)]


def setup(Solver, nx, ny, dx, dy, dtmin):
    """Build a complete Noh problem for one solver. Every object is returned and
    kept alive by the caller: the physics package stores non-owning pointers."""
    grid = Grid2d(nx, ny, dx, dy)
    nodeList = NodeList(grid.size())
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = Solver(nodeList, constants, eos, grid)
    boundary = OutflowGridBoundary2d(grid=grid)
    hydro.addBoundary(boundary)
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    x0 = 0.5 * nx * dx
    density = nodeList.getFieldDouble("density")
    energy = nodeList.getFieldDouble("specificInternalEnergy")
    velocity = nodeList.getFieldVector2d("velocity")
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            x = (i + 0.5) * dx
            vx = V0 if x < x0 else -V0          # both streams move toward x0
            density.setValue(idx, RHO0)
            energy.setValue(idx, P0 / ((GAMMA - 1) * RHO0))
            velocity.setValue(idx, Vector2d(vx, 0.0))

    return dict(grid=grid, nodeList=nodeList, constants=constants, eos=eos,
                hydro=hydro, boundary=boundary, integrator=integrator, x0=x0)


if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       animateSolver = "HLLE",
                                       cycles = 20000,
                                       nx = 200,
                                       ny = 10,
                                       dx = 0.005,
                                       dy = 0.005,
                                       tstop = 0.3,
                                       dtmin = 1e-7,
                                       plot = True)

    if animate:
        byName = dict(SOLVERS)
        assert animateSolver in byName, f"--animateSolver must be one of {list(byName)}"
        problem = setup(byName[animateSolver], nx, ny, dx, dy, dtmin)
        hydro = problem["hydro"]
        controller = Controller(integrator=problem["integrator"], periodicWork=[],
                                statStep=200, tstop=tstop)
        title = MakeTitle(controller, "time", "time")
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d, stepper=controller.Step,
                                                title=title, fieldName="density")
        AnimateGrid2d((nx, ny), update_method, extremis=[0, 5], frames=cycles, cmap="plasma",
                      lineout_axis='y', lineout_index=ny // 2)
    else:
        noh = NohSolution(GAMMA, GEOMETRY, V0, RHO0)
        rho_shock = noh.rho2
        results = []
        failures = []
        for label, Solver in SOLVERS:
            problem = setup(Solver, nx, ny, dx, dy, dtmin)
            grid, nodeList, x0 = problem["grid"], problem["nodeList"], problem["x0"]
            controller = Controller(integrator=problem["integrator"], periodicWork=[],
                                    statStep=200, tstop=tstop)
            controller.Step(cycles)
            t = controller.time
            xs = noh.shockRadius(t)                            # shock half-extent

            density = nodeList.getFieldDouble("density")
            pressure = nodeList.getFieldDouble("pressure")
            velocity = nodeList.getFieldVector2d("velocity")
            jmid = ny // 2

            # (1) exact-solution L1 density error, excluding the wall-heating band
            #     (a few cells around x0) where every scheme is anomalous.
            band = 4 * dx
            err = 0.0; n = 0
            for i in range(nx):
                idx = grid.index(i, jmid, 0)
                x = (i + 0.5) * dx
                if abs(x - x0) < band:
                    continue
                rho_e = rho_shock if abs(x - x0) < xs else RHO0
                err += abs(density[idx] - rho_e)
                n += 1
            err /= n

            # (2) plateau density and post-shock velocity, sampled between the
            #     wall-heating band and the shock front.
            plateau, vplateau = [], []
            for i in range(nx):
                idx = grid.index(i, jmid, 0)
                x = (i + 0.5) * dx
                if band < abs(x - x0) < 0.7 * xs:
                    plateau.append(density[idx])
                    vplateau.append(velocity[idx].x)
            plateau_rho = float(np.median(plateau))
            plateau_v = float(np.max(np.abs(vplateau)))

            # (3) shock-front position: outermost cell above the midpoint density.
            thresh = 0.5 * (RHO0 + rho_shock)
            x_front = 0.0
            for i in range(nx):
                idx = grid.index(i, jmid, 0)
                x = (i + 0.5) * dx
                if x > x0 and density[idx] > thresh:
                    x_front = max(x_front, x - x0)

            print(f"{label} Noh at t={t:.4f}: L1(rho)={err:.4e}  plateau rho={plateau_rho:.3f} "
                  f"(exact {rho_shock:.3f})  |v|={plateau_v:.3e}  front={x_front:.4f} "
                  f"(exact {xs:.4f})")
            if not err < 5.0e-2:
                failures.append(f"{label}: density L1 error too large: {err:.4e}")
            if not abs(plateau_rho - rho_shock) < 0.1 * rho_shock:
                failures.append(f"{label}: plateau density off: {plateau_rho:.3f} "
                                f"vs {rho_shock:.3f}")
            if not plateau_v < 5.0e-2:
                failures.append(f"{label}: post-shock velocity not stagnant: {plateau_v:.3e}")
            if not abs(x_front - xs) < 5 * dx:
                failures.append(f"{label}: shock front misplaced: {x_front:.4f} vs {xs:.4f}")

            # Fold both streams onto a single radius-from-stagnation-point axis
            # -- the outward-projected (signed) velocity component so both sides
            # show the same "negative = inbound" convention as the analytic profile.
            r_sim = np.array([abs((i + 0.5) * dx - x0) for i in range(nx)])
            rho_sim = np.array([density[grid.index(i, jmid, 0)] for i in range(nx)])
            p_sim = np.array([pressure[grid.index(i, jmid, 0)] for i in range(nx)])
            u_sim = np.array([velocity[grid.index(i, jmid, 0)].x * (1.0 if (i + 0.5) * dx > x0 else -1.0)
                              for i in range(nx)])
            results.append((label, t, r_sim, {"density": rho_sim,
                                              "pressure": p_sim,
                                              "radial velocity": u_sim}))

        assert not failures, "\n".join(failures)
        print("[ok] Noh problem matches the exact solution (density jump, stagnation, "
              "shock speed) for KT, HLLC, and HLLE.")

        if plot:
            # Each solver overshoots tstop by at most one dt, so their achieved
            # times differ negligibly; the analytic curve uses the first one's.
            label0, t_plot, r0, _ = results[0]
            r_fine = np.linspace(0.0, r0.max(), 400)
            rho_fine, u_fine, p_fine = noh.profile(t_plot, r_fine)

            fig, axes = plt.subplots(1, 3, figsize=(14, 4))
            for ax, fine, name in ((axes[0], rho_fine, "density"),
                                   (axes[1], p_fine, "pressure"),
                                   (axes[2], u_fine, "radial velocity")):
                ax.plot(r_fine, fine, "k-", lw=1, label="Noh analytic")
                for label, t, r_sim, fields in results:
                    ax.plot(r_sim, fields[name], "o", ms=2.5, alpha=0.7, label=label)
                ax.set_xlabel("r = |x-x0|"); ax.set_ylabel(name); ax.legend()
            fig.suptitle(f"Planar Noh problem at t={t_plot:.4f}: KT vs HLLC vs HLLE")
            fig.tight_layout()
            plt.show()
