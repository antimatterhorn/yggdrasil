from yggdrasil import *
import numpy as np
from Animation import *
from Physics import GridHydroHLLE2d
from Mesh import Grid2d, Geometry
from EOS import IdealGasEOS
from Boundaries import OutflowGridBoundary2d

# Cylindrical Noh problem in axisymmetric (r,z) coordinates, verified against its
# exact solution. Cold gas (rho=1, p~0) falls radially inward at speed v0,
# stagnates on the r=0 axis, and drives a shock back outward. For gamma=5/3,
# v0=1 the exact solution is:
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

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = False,
                                       cycles = 20000,
                                       nr = 200,
                                       nz = 10,
                                       dr = 0.005,
                                       dz = 0.005,
                                       tstop = 0.6,
                                       dtmin = 1e-7)

    myGrid = Grid2d(nr, nz, dr, dz, Geometry.CylindricalRZ)
    myNodeList = NodeList(nr * nz)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroHLLE2d(myNodeList, constants, eos, myGrid)

    # The r=0 axis boundary is installed automatically for RZ grids.
    box = OutflowGridBoundary2d(grid=myGrid)
    box.setFaces(["right", "top", "bottom"])
    hydro.addBoundary(box)

    integrator = RungeKutta4Integrator2d([hydro], dtmin=dtmin, verbose=False)

    density = myNodeList.getFieldDouble("density")
    energy = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")
    for j in range(nz):
        for i in range(nr):
            idx = myGrid.index(i, j, 0)
            density.setValue(idx, RHO0)
            energy.setValue(idx, P0 / ((GAMMA - 1) * RHO0))
            velocity.setValue(idx, Vector2d(-V0, 0.0))     # inward, toward r=0

    controller = Controller(integrator=integrator, periodicWork=[], statStep=200, tstop=tstop)

    if animate:
        title = MakeTitle(controller, "time", "time")
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d, stepper=controller.Step,
                                                title=title, fieldName="density")
        AnimateGrid2d((nr, nz), update_method, extremis=[0, 18], frames=cycles, cmap="plasma",
                      lineout_axis='y', lineout_index=nz // 2)
    else:
        controller.Step(cycles)
        t = controller.time

        # Exact solution.
        rho_shock = RHO0 * ((GAMMA + 1) / (GAMMA - 1)) ** 2    # = 16
        D = V0 * (GAMMA - 1) / 2.0                             # shock speed = 1/3
        rs = D * t                                             # shock radius

        def rho_exact(r):
            return rho_shock if r < rs else RHO0 * (1.0 + V0 * t / r)

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
            idx = myGrid.index(i, jmid, 0)
            r = myGrid.cellRadius(idx)
            if r < band or r > r_outer:
                continue
            err += abs(density[idx] - rho_exact(r))
            n += 1
        err /= n

        # (2) shocked plateau density and stagnation, sampled between the
        #     wall-heating band and the shock front.
        plateau, vplateau = [], []
        for i in range(nr):
            idx = myGrid.index(i, jmid, 0)
            r = myGrid.cellRadius(idx)
            if band < r < 0.7 * rs:
                plateau.append(density[idx])
                vplateau.append(velocity[idx].x)
        plateau_rho = float(np.median(plateau))
        plateau_v = float(np.max(np.abs(vplateau)))

        # (3) shock-front radius: outermost cell above the midpoint density.
        thresh = 0.5 * (rho_exact(rs * 1.01) + rho_shock)
        r_front = 0.0
        for i in range(nr):
            idx = myGrid.index(i, jmid, 0)
            if density[idx] > thresh:
                r_front = max(r_front, myGrid.cellRadius(idx))

        # (4) pre-shock convergent profile rho = rho0*(1 + v0 t / r). This is the
        #     signature of cylindrical convergence -- the planar Noh leaves the
        #     pre-shock density at rho0.
        preshock_err = 0.0
        for i in range(nr):
            idx = myGrid.index(i, jmid, 0)
            r = myGrid.cellRadius(idx)
            if 1.5 * rs < r < r_outer:
                rho_e = RHO0 * (1.0 + V0 * t / r)
                preshock_err = max(preshock_err, abs(density[idx] - rho_e) / rho_e)

        print(f"NohRZ at t={t:.4f}: L1(rho)={err:.4e}  plateau rho={plateau_rho:.3f} "
              f"(exact {rho_shock:.3f})  |v_r|={plateau_v:.3e}  front={r_front:.4f} "
              f"(exact {rs:.4f})  preshock max rel err={preshock_err:.3e}")
        assert err < 1.0, f"density L1 error too large: {err:.4e}"
        assert abs(plateau_rho - rho_shock) < 0.15 * rho_shock, \
            f"plateau density off: {plateau_rho:.3f} vs {rho_shock:.3f}"
        assert plateau_v < 0.1, f"post-shock velocity not stagnant: {plateau_v:.3e}"
        assert abs(r_front - rs) < 8 * dr, f"shock front misplaced: {r_front:.4f} vs {rs:.4f}"
        assert preshock_err < 0.02, \
            f"pre-shock convergent profile off: max rel err {preshock_err:.3e}"
        print("[ok] Cylindrical Noh matches the exact solution "
              "(r^2 compression, convergent pre-shock profile, stagnation, shock speed).")
