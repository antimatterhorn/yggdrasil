from yggdrasil import *
import sys
import math
from Physics import DEM2d, DEMConstantForce2d

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

_results = []

def check(label, condition):
    status = "PASS" if condition else "FAIL"
    print(f"    [{status}] {label}")
    _results.append(condition)


def make_two_sphere(kn, kt, gammaN, gammaT, mu, v0=1.0, offset_y=0.0):
    """Return (nodes, dem, constants, controller).

    All four must be kept alive by the caller: the integrator holds a raw
    C++ pointer to dem, and dem holds a raw reference to constants.  If
    Python GC's either before the controller finishes stepping, the process
    will segfault (same caveat as in boundary_test.py::make_hydro).
    """
    constants = MKS()
    nodes = NodeList(2)

    dem = DEM2d(nodes, constants, kn, kt, gammaN, gammaT, mu)
    integrator = RungeKutta2Integrator2d(packages=[dem], dtmin=1e-8, verbose=False)

    radius = 0.5
    m      = 1.0
    for i in range(2):
        nodes.getFieldDouble("radius").setValue(i, radius)
        nodes.getFieldDouble("mass").setValue(i, m)

    gap = 0.05
    nodes.getFieldVector2d("position").setValue(0, Vector2d(-radius - gap / 2,  offset_y / 2))
    nodes.getFieldVector2d("position").setValue(1, Vector2d( radius + gap / 2, -offset_y / 2))
    nodes.getFieldVector2d("velocity").setValue(0, Vector2d( v0, 0.0))
    nodes.getFieldVector2d("velocity").setValue(1, Vector2d(-v0, 0.0))

    controller = Controller(integrator=integrator, periodicWork=[], statStep=1000000, tstop=1.0)
    return nodes, dem, constants, controller


# ─────────────────────────────────────────────────────────────────────────────
# Test 1 — elastic head-on collision
#   No damping, no friction.  Equal masses ⟹ velocities exactly swap.
#   Momentum must be conserved to floating-point precision (Newton's 3rd law
#   is applied pair-wise in DEM).  Kinetic energy should be conserved; a small
#   numerical error from discrete contact detection is expected but < 5 %.
# ─────────────────────────────────────────────────────────────────────────────

def test_elastic_collision():
    print("Test 1: Elastic head-on collision (gammaN=0, mu=0)")
    kn, kt = 1e4, 7e3
    v0, m  = 1.0, 1.0

    nodes, _dem, _constants, controller = make_two_sphere(kn=kn, kt=kt, gammaN=0.0, gammaT=0.0, mu=0.0, v0=v0)
    controller.Step(1000000)

    vel = nodes.getFieldVector2d("velocity")
    v0f, v1f = vel[0], vel[1]

    px     = m * (v0f.x + v1f.x)
    py     = m * (v0f.y + v1f.y)
    ke0    = m * v0 ** 2                          # 2 × (½mv²)
    ke_fin = 0.5 * m * (v0f.mag2 + v1f.mag2)
    ke_err = abs(ke_fin - ke0) / ke0

    check("total px conserved  |px| < 1e-10",  abs(px) < 1e-10)
    check("total py conserved  |py| < 1e-10",  abs(py) < 1e-10)
    check("KE conserved to < 5 %",              ke_err < 0.05)
    check("particle 0 reversed  v0.x < 0",      v0f.x < 0)
    check("particle 1 reversed  v1.x > 0",      v1f.x > 0)
    print(f"      KE error = {ke_err:.4f},  px = {px:.2e},  py = {py:.2e}")
    print(f"      v0 = ({v0f.x:.4f}, {v0f.y:.4f}),  v1 = ({v1f.x:.4f}, {v1f.y:.4f})")


# ─────────────────────────────────────────────────────────────────────────────
# Test 2 — damped head-on collision, coefficient of restitution
#   For a linear spring-dashpot the analytic COR is
#       e = exp(-ζπ / √(1−ζ²))   where  ζ = γN / (2√(kn · m_red))
#   The numeric COR is |v_final| / v0 (equal masses → symmetric bounce).
# ─────────────────────────────────────────────────────────────────────────────

def test_damped_collision():
    print("Test 2: Damped head-on collision (gammaN=20, mu=0)")
    kn, kt, gammaN = 1e4, 7e3, 20.0
    v0, m = 1.0, 1.0

    m_red  = m / 2.0
    zeta   = gammaN / (2.0 * math.sqrt(kn * m_red))
    e_ana  = math.exp(-zeta * math.pi / math.sqrt(max(1.0 - zeta ** 2, 1e-12)))

    nodes, _dem, _constants, controller = make_two_sphere(kn=kn, kt=kt, gammaN=gammaN, gammaT=0.0, mu=0.0, v0=v0)
    controller.Step(1000000)

    vel = nodes.getFieldVector2d("velocity")
    v0f, v1f = vel[0], vel[1]

    px    = m * (v0f.x + v1f.x)
    ke0   = m * v0 ** 2
    ke_fin = 0.5 * m * (v0f.mag2 + v1f.mag2)
    e_num = abs(v0f.x) / v0

    check("total px conserved  |px| < 1e-10",                    abs(px) < 1e-10)
    check("energy dissipated   KE_final < KE_initial",           ke_fin < ke0)
    check(f"COR within 10 % of analytic ({e_ana:.3f})",          abs(e_num - e_ana) / e_ana < 0.10)
    print(f"      e_analytic = {e_ana:.4f},  e_numeric = {e_num:.4f}  "
          f"(error = {abs(e_num - e_ana) / e_ana:.3f})")


# ─────────────────────────────────────────────────────────────────────────────
# Test 3 — frictional oblique collision, spin generation
#   Particles are offset slightly in y so the contact is not purely normal.
#   With mu > 0 the tangential Cundall-Strack spring must produce non-zero
#   torques, leaving both particles with non-zero angular velocity.
# ─────────────────────────────────────────────────────────────────────────────

def test_frictional_spin():
    print("Test 3: Frictional oblique collision (mu=0.3) — spin generation")
    kn, kt = 1e4, 7e3
    v0     = 1.0

    nodes, _dem, _constants, controller = make_two_sphere(kn=kn, kt=kt, gammaN=5.0, gammaT=5.0,
                                                          mu=0.3, v0=v0, offset_y=0.3)
    controller.Step(1000000)

    omega = nodes.getFieldDouble("angularVelocity")
    w0, w1 = omega[0], omega[1]

    check("|ω₀| > 1e-6 after frictional contact", abs(w0) > 1e-6)
    check("|ω₁| > 1e-6 after frictional contact", abs(w1) > 1e-6)
    print(f"      ω₀ = {w0:.6f},  ω₁ = {w1:.6f}")


# ─────────────────────────────────────────────────────────────────────────────
# Optional visual — granular pile under gravity
#   DEMConstantForce runs first so the gravitational velocity kick is applied
#   before DEM reads the state.  Floor particles (bottom row) have very large
#   mass so gravity accelerates them negligibly while DEM contact forces keep
#   the pile from passing through.
# ─────────────────────────────────────────────────────────────────────────────

def run_pile(nx, ny):
    from HCPNodeGenerator import HCPNodeGenerator2d
    from Animation import AnimateScatter

    gen        = HCPNodeGenerator2d(nx, ny)
    fallers    = len(gen.positions)
    n_floor    = nx + 1
    numNodes   = fallers + n_floor

    constants  = MKS()
    nodes      = NodeList(numNodes)

    kn, kt     = 5e4, 3.5e4
    g_vec      = Vector2d(0.0, -9.8)

    body_force = DEMConstantForce2d(nodes, constants, g_vec)
    dem        = DEM2d(nodes, constants, kn, kt, 80.0, 30.0, 0.3)

    integrator = RungeKutta2Integrator2d(packages=[body_force, dem],
                                         dtmin=1e-8, verbose=False)

    radius     = 0.45
    rad_fld    = nodes.getFieldDouble("radius")
    mass_fld   = nodes.getFieldDouble("mass")
    pos_fld    = nodes.getFieldVector2d("position")

    # Falling particles from HCP generator, raised above y = 0
    for i in range(fallers):
        rad_fld.setValue(i, radius)
        mass_fld.setValue(i, 0.1)
        x = gen.positions[i][0] * (nx * radius * 1.9)
        y = gen.positions[i][1] * (ny * radius * 2.1) + 2.5
        pos_fld.setValue(i, Vector2d(x, y))

    # Floor row: large mass so body force and contact forces barely move them
    floor_mass = 1e8
    for k in range(n_floor):
        idx = fallers + k
        rad_fld.setValue(idx, radius)
        mass_fld.setValue(idx, floor_mass)
        x = (k - n_floor // 2) * radius * 2.0
        pos_fld.setValue(idx, Vector2d(x, -radius))

    controller = Controller(integrator=integrator, periodicWork=[], statStep=100)

    bounds = (-nx * radius * 1.2, nx * radius * 1.2,
              -radius * 3,        ny * radius * 2.5 + 3)
    AnimateScatter(bounds, stepper=controller, positions=pos_fld,
                   frames=400, interval=20)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate=False, nx=6, ny=6)

    test_elastic_collision()
    test_damped_collision()
    test_frictional_spin()

    n_pass  = sum(_results)
    n_total = len(_results)
    print()
    if all(_results):
        print(f"ALL {n_total} CHECKS PASSED")
    else:
        print(f"{n_total - n_pass} / {n_total} CHECKS FAILED")

    if animate:
        run_pile(nx, ny)

    sys.exit(0 if all(_results) else 1)
