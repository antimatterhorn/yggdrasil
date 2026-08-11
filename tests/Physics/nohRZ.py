from yggdrasil import *
import numpy as np
import matplotlib.pyplot as plt
from Animation import *
from Physics import GridHydroKT2d, GridHydroHLLC2d, GridHydroHLLE2d
from Mesh import Grid2d, Geometry
from EOS import IdealGasEOS
from Boundaries import OutflowGridBoundary2d
from AnalyticSolutions import NohSolution

# Cylindrical Noh problem in axisymmetric (r,z) coordinates, verified against its
# exact solution (AnalyticSolutions.py) for all three grid solvers (KT, HLLC,
# HLLE). Cold gas (rho=1, p~0) falls radially inward at speed v0, stagnates on
# the r=0 axis, and drives a shock back outward. For gamma=5/3, v0=1 the exact
# solution is:
#   * shocked region r < D t : rho = rho0*((g+1)/(g-1))^2 = 16, v = 0,
#                              p = (g-1)*rho2*v0^2/2 = 16/3
#   * unshocked region       : rho = rho0*(1 + v0 t / r), v = -v0, p ~ 0
# with shock speed D = v0*(g-1)/2 = 1/3.
#
# Two things differ from the planar Noh (noh.py): the compression is the *square*
# of the planar jump (16 vs 4) because the flow converges in two dimensions, and
# the pre-shock density is no longer uniform -- it grows as 1 + v0 t / r as mass
# piles up toward the axis. The axis suffers the usual Noh wall heating (worse in
# convergent geometry), so the plateau is checked away from r=0.

GAMMA = 5.0 / 3.0
V0 = 1.0
RHO0 = 1.0
P0 = 1.0e-6
GEOMETRY = 2                                # cylindrical: converges in r only, z is free

SOLVERS = [("KT", GridHydroKT2d),
           ("HLLC", GridHydroHLLC2d),
           ("HLLE", GridHydroHLLE2d)]


def setup(Solver, nr, nz, dr, dz, dtmin):
    """Build a complete RZ Noh problem for one solver. Every object is returned and
    kept alive by the caller: the physics package stores non-owning pointers."""
    grid = Grid2d(nr, nz, dr, dz, Geometry.CylindricalRZ)
    nodeList = NodeList(nr * nz)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = Solver(nodeList, constants, eos, grid)

    # The r=0 axis boundary is installed automatically for RZ grids.
    box = OutflowGridBoundary2d(grid=grid)
    box.setFaces(["right", "top", "bottom"])
    hydro.addBoundary(box)

    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    density = nodeList.getFieldDouble("density")
    energy = nodeList.getFieldDouble("specificInternalEnergy")
    velocity = nodeList.getFieldVector2d("velocity")
    for j in range(nz):
        for i in range(nr):
            idx = grid.index(i, j, 0)
            density.setValue(idx, RHO0)
            energy.setValue(idx, P0 / ((GAMMA - 1) * RHO0))
            velocity.setValue(idx, Vector2d(-V0, 0.0))     # inward, toward r=0

    return dict(grid=grid, nodeList=nodeList, constants=constants, eos=eos,
                hydro=hydro, boundary=box, integrator=integrator)


if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       animateSolver = "HLLE",
                                       cycles = 20000,
                                       nr = 200,
                                       nz = 10,
                                       dr = 0.005,
                                       dz = 0.005,
                                       tstop = 0.6,
                                       dtmin = 1e-7,
                                       plot = True)

    if animate:
        byName = dict(SOLVERS)
        assert animateSolver in byName, f"--animateSolver must be one of {list(byName)}"
        problem = setup(byName[animateSolver], nr, nz, dr, dz, dtmin)
        hydro = problem["hydro"]
        controller = Controller(integrator=problem["integrator"], periodicWork=[],
                                statStep=200, tstop=tstop)
        title = MakeTitle(controller, "time", "time")
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d, stepper=controller.Step,
                                                title=title, fieldName="density")
        AnimateGrid2d((nr, nz), update_method, extremis=[0, 18], frames=cycles, cmap="plasma",
                      lineout_axis='y', lineout_index=nz // 2)
    else:
        noh = NohSolution(GAMMA, GEOMETRY, V0, RHO0)
        rho_shock = noh.rho2
        results = []
        failures = []
        for label, Solver in SOLVERS:
            problem = setup(Solver, nr, nz, dr, dz, dtmin)
            grid, nodeList = problem["grid"], problem["nodeList"]
            controller = Controller(integrator=problem["integrator"], periodicWork=[],
                                    statStep=200, tstop=tstop)
            controller.Step(cycles)
            t = controller.time
            rs = noh.shockRadius(t)                                # shock radius

            def rho_exact(r):
                return rho_shock if r < rs else RHO0 * (1.0 + V0 * t / r)

            density = nodeList.getFieldDouble("density")
            pressure = nodeList.getFieldDouble("pressure")
            velocity = nodeList.getFieldVector2d("velocity")
            jmid = nz // 2
            band = 6 * dr                 # wall-heating band at the axis
            # The exact pre-shock profile assumes gas streaming in from infinity; on
            # a finite domain the zero-gradient outer boundary over-supplies, and the
            # deviation grows with radius (<1% inside r/R = 0.5, ~14% at the edge).
            # Compare only over the clean interior.
            r_outer = 0.5 * nr * dr

            # (1) exact-solution L1 density error over the clean interior.
            err = 0.0; n = 0
            for i in range(nr):
                idx = grid.index(i, jmid, 0)
                r = grid.cellRadius(idx)
                if r < band or r > r_outer:
                    continue
                err += abs(density[idx] - rho_exact(r))
                n += 1
            err /= n

            # (2) shocked plateau density and stagnation, sampled between the
            #     wall-heating band and the shock front.
            plateau, vplateau = [], []
            for i in range(nr):
                idx = grid.index(i, jmid, 0)
                r = grid.cellRadius(idx)
                if band < r < 0.7 * rs:
                    plateau.append(density[idx])
                    vplateau.append(velocity[idx].x)
            plateau_rho = float(np.median(plateau))
            plateau_v = float(np.max(np.abs(vplateau)))

            # (3) shock-front radius: outermost cell above the midpoint density.
            thresh = 0.5 * (rho_exact(rs * 1.01) + rho_shock)
            r_front = 0.0
            for i in range(nr):
                idx = grid.index(i, jmid, 0)
                if density[idx] > thresh:
                    r_front = max(r_front, grid.cellRadius(idx))

            # (4) pre-shock convergent profile rho = rho0*(1 + v0 t / r). This is the
            #     signature of cylindrical convergence -- the planar Noh leaves the
            #     pre-shock density at rho0.
            preshock_err = 0.0
            for i in range(nr):
                idx = grid.index(i, jmid, 0)
                r = grid.cellRadius(idx)
                if 1.5 * rs < r < r_outer:
                    rho_e = RHO0 * (1.0 + V0 * t / r)
                    preshock_err = max(preshock_err, abs(density[idx] - rho_e) / rho_e)

            print(f"{label} NohRZ at t={t:.4f}: L1(rho)={err:.4e}  plateau rho={plateau_rho:.3f} "
                  f"(exact {rho_shock:.3f})  |v_r|={plateau_v:.3e}  front={r_front:.4f} "
                  f"(exact {rs:.4f})  preshock max rel err={preshock_err:.3e}")
            if not err < 1.0:
                failures.append(f"{label}: density L1 error too large: {err:.4e}")
            if not abs(plateau_rho - rho_shock) < 0.15 * rho_shock:
                failures.append(f"{label}: plateau density off: {plateau_rho:.3f} "
                                f"vs {rho_shock:.3f}")
            if not plateau_v < 0.1:
                failures.append(f"{label}: post-shock velocity not stagnant: {plateau_v:.3e}")
            if not abs(r_front - rs) < 8 * dr:
                failures.append(f"{label}: shock front misplaced: {r_front:.4f} vs {rs:.4f}")
            if not preshock_err < 0.02:
                failures.append(f"{label}: pre-shock convergent profile off: "
                                f"max rel err {preshock_err:.3e}")

            r_sim = np.array([grid.cellRadius(grid.index(i, jmid, 0)) for i in range(nr)])
            rho_sim = np.array([density[grid.index(i, jmid, 0)] for i in range(nr)])
            p_sim = np.array([pressure[grid.index(i, jmid, 0)] for i in range(nr)])
            u_sim = np.array([velocity[grid.index(i, jmid, 0)].x for i in range(nr)])
            results.append((label, t, r_sim, {"density": rho_sim,
                                              "pressure": p_sim,
                                              "radial velocity": u_sim}))

        assert not failures, "\n".join(failures)
        print("[ok] Cylindrical Noh matches the exact solution "
              "(r^2 compression, convergent pre-shock profile, stagnation, shock speed) "
              "for KT, HLLC, and HLLE.")

        if plot:
            # Analytic pre-shock profile diverges as r->0 (only meaningful ahead
            # of the shock); evaluate on the simulation's own radii directly.
            # Each solver overshoots tstop by at most one dt, so their achieved
            # times differ negligibly; the analytic curve uses the first one's.
            label0, t_plot, r0, _ = results[0]
            rho_fine, u_fine, p_fine = noh.profile(t_plot, r0)

            fig, axes = plt.subplots(1, 3, figsize=(14, 4))
            for ax, fine, name in ((axes[0], rho_fine, "density"),
                                   (axes[1], p_fine, "pressure"),
                                   (axes[2], u_fine, "radial velocity")):
                ax.plot(r0, fine, "k-", lw=1, label="Noh analytic")
                for label, t, r_sim, fields in results:
                    ax.plot(r_sim, fields[name], "o", ms=2.5, alpha=0.7, label=label)
                ax.set_xlabel("r"); ax.set_ylabel(name); ax.legend()
            fig.suptitle(f"Cylindrical (r,z) Noh problem at t={t_plot:.4f}: KT vs HLLC vs HLLE")
            fig.tight_layout()
            plt.show()
