from yggdrasil import *
from Physics import GridHydroKT2d
from Mesh import Grid2d
from EOS import IdealGasEOS
from Boundaries import ReflectingGridBoundary2d

# Headless Cartesian regression: run a fixed Sod problem for a fixed number of
# cycles and print a high-precision checksum of the evolved fields. Used to
# prove the finite-volume rewrite of GridHydroBase is numerically identical to
# the previous (flux_L - flux_R)/dx formulation for Cartesian geometry. Uses the
# KT solver (HLLC/HLLE hit a pre-existing sonic-point instability at this cycle
# count regardless of the divergence formulation). Two spacings are tested: unit
# (where the rewrite is bit-identical) and non-unit (algebraically identical,
# equal to round-off).

def run(nx, ny, dx, dy, cycles=200):
    grid = Grid2d(nx, ny, dx, dy)
    nodes = NodeList(nx * ny)
    constants = MKS()
    eos = IdealGasEOS(1.4, constants)
    hydro = GridHydroKT2d(nodes, constants, eos, grid)
    hydro.addBoundary(ReflectingGridBoundary2d(grid=grid))
    integrator = RungeKutta4Integrator2d([hydro], dtmin=0.001, verbose=False)

    density = nodes.getFieldDouble("density")
    energy = nodes.getFieldDouble("specificInternalEnergy")
    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            if i < nx // 2:
                density.setValue(idx, 1.0);   energy.setValue(idx, 2.5)
            else:
                density.setValue(idx, 0.125); energy.setValue(idx, 2.0)

    Controller(integrator=integrator, periodicWork=[], statStep=100000).Step(cycles)

    velocity = nodes.getFieldVector2d("velocity")
    srho = ssie = svel = 0.0
    for i in range(nx * ny):
        srho += density[i]
        ssie += energy[i]
        svel += abs(velocity[i].x) + abs(velocity[i].y)
    return srho, ssie, svel


if __name__ == "__main__":
    for label, (nx, ny, dx, dy) in (("unit   ", (80, 16, 1.0, 1.0)),
                                     ("nonunit", (80, 16, 0.02, 0.05))):
        srho, ssie, svel = run(nx, ny, dx, dy)
        print(f"{label} sum(rho)={srho:.15e} sum(sie)={ssie:.15e} sum(|v|)={svel:.15e}")
