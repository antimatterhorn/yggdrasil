from yggdrasil import *
import numpy as np
from Animation import *
from Physics import GridHydroHLLE2d
from Mesh import Grid2d
from EOS import IdealGasEOS
from Boundaries import OutflowGridBoundary2d

# Planar Noh problem, verified against its exact solution. Two cold streams
# (rho=1, p~0) collide head-on at x0 with speed v0, driving a pair of shocks
# outward from the stagnation point. For gamma=5/3, v0=1 the exact solution is:
#   * shocked region |x-x0| < D t : rho = rho0*(g+1)/(g-1) = 4, v = 0,
#                                    p = rho0*v0^2*(g+1)/2 = 4/3
#   * unshocked region             : rho = rho0, v = +/- v0 (toward x0), p ~ 0
# with shock speed D = v0*(g-1)/2 = 1/3. The stagnation point suffers the
# well-known Noh "wall heating" over a few cells, so the plateau is checked away
# from x0.

GAMMA = 5.0 / 3.0
V0 = 1.0
RHO0 = 1.0
P0 = 1.0e-6

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       cycles = 20000,
                                       nx = 200,
                                       ny = 10,
                                       dx = 0.005,
                                       dy = 0.005,
                                       tstop = 0.3,
                                       dtmin = 1e-7)

    myGrid = Grid2d(nx, ny, dx, dy)
    myNodeList = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroHLLE2d(myNodeList, constants, eos, myGrid)
    hydro.addBoundary(OutflowGridBoundary2d(grid=myGrid))
    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    x0 = 0.5 * nx * dx
    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")
    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i, j, 0)
            x = (i + 0.5) * dx
            vx = V0 if x < x0 else -V0          # both streams move toward x0
            density.setValue(idx, RHO0)
            energy.setValue(idx, P0 / ((GAMMA - 1) * RHO0))
            velocity.setValue(idx, Vector2d(vx, 0.0))

    controller = Controller(integrator=integrator, periodicWork=[], statStep=200, tstop=tstop)

    if animate:
        title = MakeTitle(controller, "time", "time")
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d, stepper=controller.Step,
                                                title=title, fieldName="density")
        AnimateGrid2d((nx, ny), update_method, extremis=[0, 5], frames=cycles, cmap="plasma")
    else:
        controller.Step(cycles)
        t = controller.time

        # Exact solution constants.
        rho_shock = RHO0 * (GAMMA + 1) / (GAMMA - 1)      # = 4
        p_shock = RHO0 * V0 * V0 * (GAMMA + 1) / 2.0       # = 4/3
        D = V0 * (GAMMA - 1) / 2.0                         # shock speed = 1/3
        xs = D * t                                         # shock half-extent

        pressure = myNodeList.getFieldDouble("pressure")
        jmid = ny // 2

        # (1) exact-solution L1 density error, excluding the wall-heating band
        #     (a few cells around x0) where every scheme is anomalous.
        band = 4 * dx
        err = 0.0; n = 0
        for i in range(nx):
            idx = myGrid.index(i, jmid, 0)
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
            idx = myGrid.index(i, jmid, 0)
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
            idx = myGrid.index(i, jmid, 0)
            x = (i + 0.5) * dx
            if x > x0 and density[idx] > thresh:
                x_front = max(x_front, x - x0)

        print(f"Noh at t={t:.4f}: L1(rho)={err:.4e}  plateau rho={plateau_rho:.3f} "
              f"(exact {rho_shock:.3f})  |v|={plateau_v:.3e}  front={x_front:.4f} "
              f"(exact {xs:.4f})")
        assert err < 5.0e-2, f"density L1 error too large: {err:.4e}"
        assert abs(plateau_rho - rho_shock) < 0.1 * rho_shock, "plateau density off"
        assert plateau_v < 5.0e-2, f"post-shock velocity not stagnant: {plateau_v:.3e}"
        assert abs(x_front - xs) < 5 * dx, f"shock front misplaced: {x_front:.4f} vs {xs:.4f}"
        print("[ok] Noh problem matches the exact solution (density jump, stagnation, shock speed).")
