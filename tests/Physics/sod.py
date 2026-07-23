from yggdrasil import *
import numpy as np
import matplotlib.pyplot as plt
from Animation import *
from Physics import GridHydroKT2d
from Mesh import Grid2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d
from AnalyticSolutions import SodSolution

# Sod shock tube, verified against the exact Riemann solution (AnalyticSolutions.py).
# The classic initial state (rho,u,p) = (1,0,1) | (0.125,0,0.1) with gamma=1.4
# develops a left rarefaction, a contact, and a right shock; the exact
# self-similar solution is sampled at the achieved time and compared to the
# numerical profile (L1 error), and plotted alongside it when plot=True.

GAMMA = 1.4

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       cycles = 20000,
                                       nx = 200,
                                       ny = 10,
                                       dx = 0.005,
                                       dy = 0.005,
                                       tstop = 0.2,
                                       dtmin = 1e-7,
                                       plot = True)

    myGrid = Grid2d(nx, ny, dx, dy)
    myNodeList = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroKT2d(myNodeList, constants, eos, myGrid)
    hydro.addBoundary(ReflectingGridBoundary2d(grid=myGrid))
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    # Classic Sod state: (rho,p) = (1,1) left of x0, (0.125,0.1) right.
    x0 = 0.5 * nx * dx
    rhoL, pL, rhoR, pR = 1.0, 1.0, 0.125, 0.1
    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            if (i + 0.5) * dx < x0:
                density.setValue(idx, rhoL); energy.setValue(idx, pL / ((GAMMA - 1) * rhoL))
            else:
                density.setValue(idx, rhoR); energy.setValue(idx, pR / ((GAMMA - 1) * rhoR))

    controller = Controller(integrator=integrator, periodicWork=[], statStep=200, tstop=tstop)

    if animate:
        title = MakeTitle(controller, "time", "time")
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d, stepper=controller.Step,
                                                title=title, fieldName="density")
        AnimateGrid2d((nx, ny), update_method, extremis=[0, 1.1], frames=cycles, cmap="plasma")
    else:
        controller.Step(cycles)
        t = controller.time
        sod = SodSolution(GAMMA, rhoL, 0.0, pL, rhoR, 0.0, pR)

        pressure = myNodeList.getFieldDouble("pressure")
        velocity = myNodeList.getFieldVector2d("velocity")
        jmid = ny // 2
        x_sim = np.array([(i + 0.5) * dx for i in range(nx)])
        rho_sim = np.array([density[myGrid.index(i, jmid, 0)] for i in range(nx)])
        p_sim = np.array([pressure[myGrid.index(i, jmid, 0)] for i in range(nx)])
        u_sim = np.array([velocity[myGrid.index(i, jmid, 0)].x for i in range(nx)])

        rho_ex, u_ex, p_ex = sod.profile(t, x_sim, x0=x0)
        err_rho = np.mean(np.abs(rho_sim - rho_ex))
        err_u = np.mean(np.abs(u_sim - u_ex))
        err_p = np.mean(np.abs(p_sim - p_ex))

        print(f"Sod at t={t:.4f}: L1 errors  rho={err_rho:.4e}  u={err_u:.4e}  p={err_p:.4e}")
        assert err_rho < 2.0e-2, f"density L1 error too large: {err_rho:.4e}"
        assert err_u < 2.0e-2, f"velocity L1 error too large: {err_u:.4e}"
        assert err_p < 2.0e-2, f"pressure L1 error too large: {err_p:.4e}"
        print("[ok] Sod shock tube matches the exact Riemann solution.")

        if plot:
            x_fine = np.linspace(x_sim.min(), x_sim.max(), 400)
            rho_fine, u_fine, p_fine = sod.profile(t, x_fine, x0=x0)

            fig, axes = plt.subplots(1, 3, figsize=(14, 4))
            for ax, sim, fine, label in ((axes[0], rho_sim, rho_fine, "density"),
                                         (axes[1], p_sim, p_fine, "pressure"),
                                         (axes[2], u_sim, u_fine, "velocity")):
                ax.plot(x_sim, sim, "o", ms=3, label="simulation")
                ax.plot(x_fine, fine, "-", label="exact Riemann")
                ax.set_xlabel("x"); ax.set_ylabel(label); ax.legend()
            fig.suptitle(f"Sod shock tube at t={t:.4f}, cycle={controller.cycle}")
            fig.tight_layout()
            plt.show()
