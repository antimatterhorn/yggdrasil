# Registry of ATS test modules. See ats.py for the analytic/snapshot contract.
tests = {
    # analytic (self-verifying against exact/closed-form solutions)
    "sod",             # Sod shock tube, KT + HLLC + HLLE vs exact Riemann
    "noh",             # planar Noh, exact stagnation-shock solution
    "sedov",           # Sedov-Taylor blast, self-similar radius law
    "sedovRZ",         # cylindrical (r,z) Sedov, per-radian mass + expansion
    "kinetics",        # free particles, constant-velocity closed form
    "kinetics_grav",   # projectile motion under constant gravity, closed form
    "orbit",           # Kepler orbit: conserved energy/momentum, period
    "kelvin-helmholtz",  # shear instability: perturbation must grow

    # snapshot (regression; no closed form)
    "plinko",          # collider maze; deterministic return to same state
    "nbody",           # chaotic 3-body
    "treeGravity",     # random N-body (Barnes-Hut)
    "wave_bounds",     # driven wave with obstacles
    "reactionDiffusion",  # rock-paper-scissors pattern formation
}
