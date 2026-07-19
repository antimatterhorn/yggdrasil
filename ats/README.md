# ATS

Automated Test System: end-to-end scripts that each build and run a full Yggdrasil simulation and check its result. This is distinct from `tests/`, which houses smaller unit-style tests — `ats/` exercises the whole simulation stack (mesh, physics package, integrator, controller) the way a real user script would.

## Two kinds of test

Each test's `run()` returns a dict whose `mode` tells the harness how to judge it:

- **`analytic`** — self-verifying. The test compares its own output against an **exact / closed-form solution** (with its own tolerances) and returns one pass/fail per check. These prove *correctness*, not merely that the output is unchanged. They need no stored reference.

  ```python
  return {"mode": "analytic", "checks": [(label, ok_bool, detail_str), ...]}
  ```

- **`snapshot`** — regression. For problems with no closed form (chaotic N-body, pattern formation, a collider maze), the returned values are diffed against `ats_reference.txt`. This only proves the output hasn't drifted. Comparison uses a **relative tolerance** (not exact equality) because OpenMP reductions aren't bit-reproducible run to run.

  ```python
  return {"mode": "snapshot", "values": [float, ...]}
  ```

A bare list return is still accepted and treated as a snapshot.

## How it works

- [yggdrasil.py](yggdrasil.py) wires `sys.path` to the built C++ modules and re-exports everything for `from yggdrasil import *`.
- [scripts.py](scripts.py) lists the test modules to run.
- [ats.py](ats.py) is the runner: it imports each test, calls `run()`, and — per `mode` — either reports the analytic checks or diffs the snapshot values against `ats_reference.txt`, with colorized pass/fail.
- [genats.py](genats.py) regenerates `ats_reference.txt`. Only **snapshot** tests are recorded; analytic tests are skipped (they carry their own reference — the exact solution).

## Contents

| Script | Mode | Verifies against |
|---|---|---|
| [sod.py](sod.py) | analytic | Exact Riemann solution — KT, HLLC, and HLLE lineouts vs the exact rarefaction/contact/shock (L1 error). Replaces the old `hllcsod`/`hllesod` snapshots. |
| [noh.py](noh.py) | analytic | Exact planar Noh: stagnation density jump `(γ+1)/(γ-1)`, `v=0` plateau, shock speed `(γ-1)/2`. |
| [sedov.py](sedov.py) | analytic | Sedov-Taylor self-similar law `R ∝ t^(1/2)` (2D → cylindrical) + energy conservation. |
| [sedovRZ.py](sedovRZ.py) | analytic | Cylindrical `(r,z)` blast: finite, expanding, per-radian mass conserved. Exercises the RZ metric, `p/r` source, and auto axis boundary. |
| [kinetics.py](kinetics.py) | analytic | Free particles: constant-velocity motion `x = x0 + v0·t` to machine precision. |
| [kinetics_grav.py](kinetics_grav.py) | analytic | Projectile motion under constant gravity `x = x0 + v0·t + ½g·t²` (RK2 integrates the linear ODE exactly). |
| [orbit.py](orbit.py) | analytic | Kepler orbit: conserved specific energy & angular momentum, and semi-major axis `a = -μ/(2E)` vs the analytical value. |
| [kelvin-helmholtz.py](kelvin-helmholtz.py) | analytic | Instability *grows*: the interface density spread climbs from 0 (flat) to a clear finite amplitude, and stays bounded. |
| [plinko.py](plinko.py) | snapshot | Particles falling through a staggered collider maze — deterministic but effectively chaotic (was the old `kinetics_grav`). |
| [nbody.py](nbody.py) | snapshot | Chaotic three-body gravity. |
| [treeGravity.py](treeGravity.py) | snapshot | Random N-body via Barnes-Hut tree gravity. |
| [wave_bounds.py](wave_bounds.py) | snapshot | Driven scalar wave with periodic + Dirichlet-obstacle boundaries. |
| [reactionDiffusion.py](reactionDiffusion.py) | snapshot | Rock-paper-scissors reaction-diffusion pattern formation. |
| [ats.py](ats.py) | — | Test runner (analytic checks + snapshot diff). |
| [genats.py](genats.py) | — | Regenerates `ats_reference.txt` (snapshot tests only). |
| [ats_reference.txt](ats_reference.txt) | — | Stored golden output for the snapshot tests. |
| [scripts.py](scripts.py) | — | Registry of test modules. |
| [yggdrasil.py](yggdrasil.py) | — | `sys.path` bootstrap + `import *` re-export. |
