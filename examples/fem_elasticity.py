from yggdrasil import *
from Mesh import FEMesh2d
from Physics import FEMLinearElasticity2d, ConstantForce2d
from Materials import IsotropicLinearElastic2d, PlaneCondition
from Boundaries import MotionConstraint2d
from Utilities import SiloDump

# Cantilever beam test using example.obj
# Left edge is pinned, gravity pulls remaining nodes down.
# Outputs silo files every dumpCycle cycles for visualization.

if __name__ == "__main__":
    commandLine = CommandLineArguments(
        cycles    = 5000,
        dtmin     = 1e-5,
        dumpCycle = 50,
        animate   = False,
    )

    # -----------------------------------------------------------------------
    # Mesh
    # -----------------------------------------------------------------------
    mesh = FEMesh2d()
    mesh.buildFromObj("example.obj", axes="(x,z)")

    nodes    = mesh.getNodes()
    numNodes = len(nodes)
    print(f"Mesh: {numNodes} nodes, {len(mesh.getElementConnectivity())} elements")

    myNodeList = NodeList(numNodes)
    constants  = MKS()

    # -----------------------------------------------------------------------
    # Material (soft, normalized)
    # E=1e4 Pa, rho=1 kg/m³  →  c_p ≈ 116 m/s
    # Min edge h_min ≈ 0.008 m → CFL limit dt_crit ≈ 4.8e-5 s (Heun stability)
    # Running at dtmin=1e-5 gives ω·dt ≈ 0.29 << √2 ≈ 1.41 — well within stability
    # -----------------------------------------------------------------------
    E   = 1e4
    nu  = 0.3
    rho = 1.0

    material = IsotropicLinearElastic2d(E, nu, PlaneCondition.Stress)

    alpha = 5.0
    beta  = 1e-5

    elasticity = FEMLinearElasticity2d(myNodeList, constants, mesh, material,
                                       rho, alpha, beta)

    comp_vec = Vector2d(-100.0, 0.0)
    compress = ConstantForce2d(myNodeList, constants, comp_vec)

    # -----------------------------------------------------------------------
    # Pin the left edge (x < x_min + 0.1) — fixed-wall cantilever support
    # -----------------------------------------------------------------------
    x_min      = min(nodes[i].x for i in range(numNodes))
    pin_thresh = x_min + 0.1
    pinned_ids = [i for i in range(numNodes) if nodes[i].x < pin_thresh]
    print(f"Pinning {len(pinned_ids)} nodes at left edge (x < {pin_thresh:.3f})")

    constraint = MotionConstraint2d(myNodeList, pinned_ids, Vector2d(1, 1))
    elasticity.addBoundary(constraint)

    # -----------------------------------------------------------------------
    # Output
    # -----------------------------------------------------------------------
    dump = SiloDump(
        baseName  = "fem_elasticity",
        nodeList  = myNodeList,
        fieldNames = ["displacement", "vonMises", "sigmaXX", "sigmaYY", "sigmaXY"],
        dumpCycle = dumpCycle,
    )

    # -----------------------------------------------------------------------
    # Integrate
    # -----------------------------------------------------------------------
    integrator = RungeKutta2Integrator2d([compress, elasticity], dtmin=dtmin, verbose=False)
    controller = Controller(integrator=integrator, periodicWork=[dump], statStep=50)

    print(f"Running {cycles} cycles...")
    controller.Step(cycles)
    print(f"Done. Final time = {controller.time:.4e} s, dt = {controller.dt:.4e} s")
