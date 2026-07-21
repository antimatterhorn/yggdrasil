from yggdrasil import *

# Free particles under zero force: each moves in a straight line at constant
# velocity, x_i(t) = x0_i + v0_i * t. Verified against that exact solution at
# the achieved time.


def run():
    numNodes = 2
    constants = MKS()
    nodes = NodeList(numNodes)

    constantForce = ConstantForce2d(nodes, constants, Vector2d(0, 0))   # zero force
    kinetics = Kinetics2d(nodes, constants)
    integrator = RungeKutta4Integrator2d(packages=[constantForce, kinetics],
                                         dtmin=0.01, verbose=False)

    pos = nodes.position
    vel = nodes.velocity
    x0 = [Vector2d(-10 + i / numNodes * 30, -0.5 + i / numNodes * 1) for i in range(numNodes)]
    v0 = [Vector2d(3 * (-1) ** i, 0) for i in range(numNodes)]
    for i in range(numNodes):
        pos[i] = x0[i]
        vel[i] = v0[i]

    Controller(integrator=integrator, periodicWork=[], statStep=100000).Step(35)
    t = integrator.time

    err = 0.0
    for i in range(numNodes):
        ex = x0[i].x + v0[i].x * t
        ey = x0[i].y + v0[i].y * t
        err = max(err, abs(pos[i].x - ex), abs(pos[i].y - ey))

    return {"mode": "analytic", "checks": [
        ("constant-velocity motion", bool(err < 1e-10),
         f"max|pos - (x0+v0 t)|={err:.2e} at t={t:.3f}"),
    ]}


if __name__ == "__main__":
    print(run())
