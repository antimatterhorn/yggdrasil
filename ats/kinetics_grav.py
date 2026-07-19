from yggdrasil import *

# Projectile motion under constant gravity, verified against the closed-form
# solution: x(t) = x0 + v0*t + 0.5*g*t^2. With a constant acceleration the ODE
# is linear, so RK2 integrates it exactly (up to round-off) -- the numerical
# trajectory must match the analytical parabola to machine precision. (The
# collider-maze "plinko" version now lives in plinko.py.)

G = -9.8


def run():
    numNodes = 5
    constants = MKS()
    nodes = NodeList(numNodes)

    gravVec = Vector2d(0, G)
    constantForce = ConstantForce2d(nodes, constants, gravVec)
    kinetics = Kinetics2d(nodes, constants)
    integrator = RungeKutta2Integrator2d(packages=[constantForce, kinetics],
                                         dtmin=0.01, verbose=False)

    mass = nodes.getFieldDouble("mass")
    pos = nodes.getFieldVector2d("position")
    vel = nodes.getFieldVector2d("velocity")
    x0, v0 = [], []
    for i in range(numNodes):
        mass[i] = 1.0
        p = Vector2d(-4 + 2 * i, 0.0)
        v = Vector2d(0.5 * i, 5.0 - i)          # varied launch velocities
        pos[i] = p; vel[i] = v
        x0.append(p); v0.append(v)

    Controller(integrator=integrator, periodicWork=[], statStep=100000).Step(120)
    t = integrator.time

    err_pos = err_vel = 0.0
    for i in range(numNodes):
        ex = x0[i].x + v0[i].x * t
        ey = x0[i].y + v0[i].y * t + 0.5 * G * t * t
        evy = v0[i].y + G * t
        err_pos = max(err_pos, abs(pos[i].x - ex), abs(pos[i].y - ey))
        err_vel = max(err_vel, abs(vel[i].x - v0[i].x), abs(vel[i].y - evy))

    return {"mode": "analytic", "checks": [
        ("position vs x0+v0 t+1/2 g t^2", bool(err_pos < 1e-9), f"max err={err_pos:.2e} at t={t:.3f}"),
        ("velocity vs v0+g t", bool(err_vel < 1e-9), f"max err={err_vel:.2e}"),
    ]}


if __name__ == "__main__":
    print(run())
