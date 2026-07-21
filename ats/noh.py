from yggdrasil import *
import numpy as np

# Planar Noh problem, verified against its exact solution. Two cold streams
# (rho=1, p~0) collide at x0 with speed v0=1, gamma=5/3, driving a pair of
# shocks outward. Exact: shocked plateau rho=rho0*(g+1)/(g-1)=4, v=0,
# p=rho0*v0^2*(g+1)/2=4/3, shock speed D=v0*(g-1)/2=1/3. The stagnation point
# suffers Noh wall-heating over a few cells, so the plateau is checked away from
# x0.

GAMMA = 5.0 / 3.0
V0, RHO0, P0 = 1.0, 1.0, 1.0e-6


def run():
    nx, ny, dx, tstop = 200, 6, 0.005, 0.3
    grid = Grid2d(nx, ny, dx, dx)
    nodes = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(GAMMA, constants)
    hydro = GridHydroHLLE2d(nodes, constants, eos, grid)
    hydro.addBoundary(OutflowGridBoundary2d(grid=grid))
    integrator = RungeKutta4Integrator2d([hydro], dtmin=1e-7, verbose=False)

    x0 = 0.5 * nx * dx
    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    velocity = nodes.getFieldVector2d("velocity")
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            density[idx] = RHO0
            energy[idx] = P0 / ((GAMMA - 1) * RHO0)
            velocity[idx] = Vector2d(V0 if (i + 0.5) * dx < x0 else -V0, 0.0)

    Controller(integrator=integrator, periodicWork=[], statStep=100000, tstop=tstop).Step(100000)
    t = integrator.time

    rho_shock = RHO0 * (GAMMA + 1) / (GAMMA - 1)      # 4
    D = V0 * (GAMMA - 1) / 2.0                         # 1/3
    xs = D * t
    band = 4 * dx
    jmid = ny // 2

    err = 0.0; n = 0
    plateau, vplateau = [], []
    x_front = 0.0
    for i in range(nx):
        idx = grid.index(i, jmid, 0)
        x = (i + 0.5) * dx
        if abs(x - x0) >= band:
            err += abs(density[idx] - (rho_shock if abs(x - x0) < xs else RHO0)); n += 1
        if band < abs(x - x0) < 0.7 * xs:
            plateau.append(density[idx]); vplateau.append(velocity[idx].x)
        if x > x0 and density[idx] > 0.5 * (RHO0 + rho_shock):
            x_front = max(x_front, x - x0)
    err /= n
    plateau_rho = float(np.median(plateau))
    plateau_v = float(np.max(np.abs(vplateau)))

    return {"mode": "analytic", "checks": [
        ("L1 density vs exact", bool(err < 5.0e-2), f"{err:.2e} < 5e-2"),
        ("shocked plateau density", bool(abs(plateau_rho - rho_shock) < 0.1 * rho_shock),
         f"{plateau_rho:.3f} vs exact {rho_shock:.3f}"),
        ("post-shock stagnation", bool(plateau_v < 5.0e-2), f"|v|={plateau_v:.2e} < 5e-2"),
        ("shock front position", bool(abs(x_front - xs) < 5 * dx),
         f"{x_front:.4f} vs exact {xs:.4f}"),
    ]}


if __name__ == "__main__":
    print(run())
