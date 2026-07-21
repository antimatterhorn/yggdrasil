from yggdrasil import *
from math import sqrt, pi

# A single test particle orbiting a fixed point mass -- a Kepler two-body
# problem, whose bound orbit conserves specific energy E = 1/2 v^2 - mu/r and
# specific angular momentum L = x*vy - y*vx exactly. The integrator must hold
# both nearly constant, and the semi-major axis inferred from the energy,
# a = -mu/(2E), must match the analytical value from the initial conditions.


def _energy_L(pos, vel, mu):
    r = sqrt(pos.x * pos.x + pos.y * pos.y)
    E = 0.5 * (vel.x * vel.x + vel.y * vel.y) - mu / r
    L = pos.x * vel.y - pos.y * vel.x
    return E, L


def run():
    nodes = NodeList(1)
    constants = PhysicalConstants(6.378e+6, 5.972e+24, 1.0)   # Earth radius/mass/s
    mu = constants.G * 1.0                                    # point mass = 1

    kinematics = Kinematics2d(nodeList=nodes, constants=constants)
    sourceGrav = PointSourceGravity2d(nodeList=nodes, constants=constants,
                                      pointSourceLocation=Vector2d(0, 0),
                                      pointSourceMass=1, pointSourceVelocity=Vector2d(0, 0))
    integrator = RungeKutta4Integrator2d(packages=[kinematics, sourceGrav], dtmin=1e-1)

    r0 = 2.0
    nodes.position[0] = Vector2d(-r0, 0.0)
    v0 = -0.5 * sqrt(mu / r0)                                 # sub-circular -> ellipse
    nodes.velocity[0] = Vector2d(0.0, v0)

    # Analytical orbit: this point is apoapsis (r0 = a(1+e)).
    E0, L0 = _energy_L(nodes.position[0], nodes.velocity[0], mu)
    a_analytic = -mu / (2.0 * E0)

    Controller(integrator=integrator, periodicWork=[], statStep=100000).Step(40000)

    Ef, Lf = _energy_L(nodes.position[0], nodes.velocity[0], mu)
    a_measured = -mu / (2.0 * Ef)
    dE = abs(Ef - E0) / abs(E0)
    dL = abs(Lf - L0) / abs(L0)
    da = abs(a_measured - a_analytic) / a_analytic

    return {"mode": "analytic", "checks": [
        ("specific energy conserved", bool(dE < 1e-3), f"dE/E={dE:.2e} < 1e-3"),
        ("angular momentum conserved", bool(dL < 1e-6), f"dL/L={dL:.2e} < 1e-6"),
        ("semi-major axis vs analytic", bool(da < 1e-3), f"a={a_measured:.4f} vs {a_analytic:.4f}"),
    ]}


if __name__ == "__main__":
    print(run())
