# ATS

Automated Test System: end-to-end regression scripts that each build and run a full Yggdrasil simulation, then compare its numerical output against a stored reference. This is distinct from `tests/`, which houses smaller unit-style tests — `ats/` exercises the whole simulation stack (mesh, physics package, integrator, controller) the way a real user script would.

## How it works

- [yggdrasil.py](yggdrasil.py) adds the built C++/pybind11 module directories to `sys.path` and re-exports every Yggdrasil submodule (`Mesh`, `Physics`, `EOS`, `Integrators`, etc.) as one `from yggdrasil import *` convenience import used by every test script.
- [scripts.py](scripts.py) lists the names of the test modules to run: `kinetics`, `kinetics_grav`, `wave_bounds`, `nbody`, `hllesod`, `hllcsod`, `orbit`, `reactionDiffusion`, `kelvin-helmholtz`.
- Each test script exposes a `run()` function that sets up and executes a simulation, then returns a flat list of numbers (times, positions, densities, etc.) summarizing the result.
- [genats.py](genats.py) imports every module listed in `scripts.py`, calls its `run()`, and writes the returned values to [ats_reference.txt](ats_reference.txt) — one `# name` header per test followed by its `repr()`'d output values. This regenerates the "golden" reference output.
- [ats.py](ats.py) is the actual test runner: it re-loads `ats_reference.txt`, re-runs every test in `scripts.py`, and diffs each script's fresh `run()` output against the stored reference, printing a colorized pass/fail summary with a unified diff for any mismatches.

## Contents

| Script | Description |
|---|---|
| [ats.py](ats.py) | Test runner — parses `ats_reference.txt`, re-runs every registered test, and diffs actual vs. expected output with colorized pass/fail reporting |
| [ats_reference.txt](ats_reference.txt) | Stored "golden" reference output (per-test lists of `repr()`'d numbers) that `ats.py` diffs fresh simulation runs against |
| [genats.py](genats.py) | Regenerates `ats_reference.txt` by running every test in `scripts.py` and recording its `run()` output |
| [hllcsod.py](hllcsod.py) | 2D Sod shock tube on a 100×20 grid using `GridHydroKT2d` (Kurganov-Tadmor flux) with an ideal-gas EOS and reflecting boundaries, tracking density along a horizontal lineout |
| [hllesod.py](hllesod.py) | Same 2D Sod shock tube setup as `hllcsod.py` but using `GridHydroHLLE2d` (HLLE Riemann solver) instead |
| [kelvin-helmholtz.py](kelvin-helmholtz.py) | Kelvin-Helmholtz shear instability on a 100×50 periodic grid with `GridHydroKT2d`, seeded with a sinusoidal velocity perturbation between two density/velocity bands, tracking a vertical density lineout |
| [kinetics.py](kinetics.py) | Minimal two-particle kinematics test using `ConstantForce2d` (zero force) and `Kinetics2d` to verify basic free-particle motion under `RungeKutta4Integrator2d` |
| [kinetics_grav.py](kinetics_grav.py) | Particles falling under constant gravity (`ConstantForce2d`, g=-9.8) through a field of `SphereCollider2d`/`BoxCollider2d` obstacles, integrated with `RungeKutta2Integrator2d` and tracking particle heights |
| [nbody.py](nbody.py) | Three-body gravitational simulation using `NBodyGravity2d` (Plummer-softened) and `Kinematics2d`, tracking the y-positions of all bodies over time |
| [orbit.py](orbit.py) | Single test-particle orbit around a fixed point mass using `PointSourceGravity2d`, initialized with a circular-orbit velocity and integrated over multiple orbital periods |
| [reactionDiffusion.py](reactionDiffusion.py) | Three-species reaction-diffusion pattern formation on a 100×100 periodic grid using the `ReactionDiffusion` physics package, seeded with a random species assignment |
| [scripts.py](scripts.py) | Registry set (`tests`) naming every test module that `ats.py` and `genats.py` should import and run |
| [wave_bounds.py](wave_bounds.py) | 2D scalar wave equation (`WaveEquation2d`) on a 100×100 grid with periodic and Dirichlet (spherical obstacle) boundaries, driven by oscillating point sources, tracking a grid of max-amplitude values |
| [yggdrasil.py](yggdrasil.py) | Shared bootstrap module that wires up `sys.path` to the built C++ extension modules and re-exports all Yggdrasil submodules for `import *` use in every other ATS script |
