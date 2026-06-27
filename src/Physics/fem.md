# FEM Structural Mechanics — Implementation Notes

## What Was Built

A complete FEM linear elasticity physics package for structural mechanics simulations. The implementation lives across three modules:

### `src/Materials/`
New module (added to root `CMakeLists.txt`). Contains:
- `material.hh` — abstract base `FEMMaterial<dim>` with pure virtual `materialMatrix()` returning an Eigen::MatrixXd (the constitutive D matrix)
- `linearElastic.cc` — `IsotropicLinearElastic<dim>`: E, ν, PlaneCondition (Stress or Strain). Builds the standard 3×3 (2D) or 6×6 (3D Voigt) D matrix.
- `Materials_PYB11.py`, `material.py`, `linearElastic.py` — Python bindings. No trampoline needed for the abstract base (it's only subclassed from C++); removing `pure_virtual` injection was required to avoid a `typename`/`int` template mismatch in the generated code.

### `src/Mesh/element.hh` — Additions
Three new methods added to `Element<dim>` (pure virtual) and implemented on `TriangleElement` and `QuadElement`:
- `computeStructuralStiffnessMatrix(positions, D)` — returns the element K_e (6×6 CST, 8×8 Q4 via 2×2 Gauss). For Q4, uses `std::abs(detJ)` everywhere — OBJ meshes often have mixed CW/CCW quad winding, and `|det| < 1e-15` skips degenerate Gauss points instead of throwing.
- `computeStrain(positions, ue)` — returns [ε_xx, ε_yy, γ_xy]. For Q4, Gauss-point averaged with the same degenerate-point skip.
- `computeLumpedMassMatrix(positions)` — already existed; patched to use `std::abs(detJ)` for mass integration so CW elements contribute positive mass.

**Why `std::abs(detJ)`:** For a CW element, `J⁻¹ = adj(J)/det(J)` has a sign flip, which propagates into B. Since `Ke = Bᵀ DB |det|`, the sign cancels in B (appears squared). For the mass, shape functions N are sign-independent; only the area measure `|det|` matters.

### `src/Physics/femLinearElasticity.cc`
The main physics package. Key design points:

**Enrolled fields:** `position`, `velocity`, `displacement` (Vector); `vonMises`, `sigmaXX`, `sigmaYY`, `sigmaXY` (scalar). State fields: `position`, `velocity`.

**`ZeroTimeInitialize`:** Copies mesh node positions into the NodeList `position` field AND into `refPositions` (member `std::vector<Vector>`). Assembles lumped nodal mass. Then calls `UpdateState()` and `InitializeBoundaries()` — the MotionConstraint's `ZeroTimeInitialize` must fire after positions are set so its `copyState` snapshot is valid.

**`EvaluateDerivatives`:** Assembles elastic forces as `f = −K_e · u_e` where `u_e = pos − refPos` (displacement from reference). Rayleigh damping: `f_damp = −(α·m + β·k_diag) · v` (diagonal approximation). Acceleration: `a = (f_elastic + f_damp)/m + bodyForce`. Derivative of position: `dxdt = v + dt·a` — this matches the framework's RK2 convention (see `constantGravity.cc` and `waveEquation.cc`; the RK2 integrator passes `dt=0` for stage 1 and `dt=actual` for stage 2).

**`FinalChecks`:** Updates `displacement` field; computes per-element strain via `computeStrain`, stress via `D·ε`, nodal-averaged `sigmaXX/YY/XY`, and von Mises `σ_vm = √(σ_xx² − σ_xx·σ_yy + σ_yy² + 3τ²)`.

**`EstimateTimestep`:** Returns `0.5 * min_i(√(m_i / k_diag_i))`. This is ~4× more conservative than the true Heun stability limit (`dt < √2 / ω_max = √2·h/c`).

**`assembleNodalMass` (private):** Iterates elements, calls `computeLumpedMassMatrix(meshPos)` scaled by `rho`, accumulates into the NodeList `mass` field. Note: `mass` is auto-enrolled by `Physics::VerifyFields`, no explicit enrollment needed.

### `src/Physics/femLinearElasticity.py` + `Physics_PYB11.py`
Python bindings follow the standard pattern. `Physics_PYB11.py` includes `material.hh`, `linearElastic.cc`, and `femLinearElasticity.cc` in its include list.

### `tests/fem_elasticity_test.py`
Cantilever test: loads `example.obj` with `axes="(x,z)"` (3687 nodes, 3374 quads), pins left-edge nodes (`x < x_min + 0.1`), applies downward gravity. Uses `SiloDump` to write fields every `dumpCycle` cycles for visualization.

## Stability / CFL Notes

The integrator's `dtmin` parameter is a **hard floor** — `VoteDt` computes `max(EstimateTimestep(), dtmin)`. If `dtmin` is larger than the actual CFL limit, the simulation will run above CFL and blow up.

For `example.obj` (min edge h ≈ 0.008 m):
- **E=1e6, ρ=1:** c_p ≈ 1160 m/s → CFL limit ≈ 4.8e-6 s → **unstable** at dtmin=1e-5 (ω·dt ≈ 2.9 > √2)
- **E=1e4, ρ=1:** c_p ≈ 116 m/s → CFL limit ≈ 4.8e-5 s → **stable** at dtmin=1e-5 (ω·dt ≈ 0.29)

The test currently uses E=1e4.

## Pending: Body Force Separation

**Current state (undesirable):** `bodyForce` is a fixed `Lin::Vector<dim>` member of `FEMLinearElasticity`, passed in the constructor and added as a constant acceleration to every node every step. This prevents time-varying or spatially varying body forces (electromagnetic, centrifugal, contact, non-uniform gravity, etc.).

**Desired direction:** Remove `bodyForce` from `FEMLinearElasticity` entirely and handle external body loads through a separate physics package (similar to how `ConstantGravity` already exists as a standalone `Kinematics` package for SPH-style simulations). For FEM this could mean:

- A `FEMBodyForce<dim>` base class (or reuse the existing `ConstantGravity`/`Kinematics` pattern) that accumulates per-node forces into a field (e.g., `"bodyForce"` Vector field on the NodeList).
- `FEMLinearElasticity::EvaluateDerivatives` reads that field instead of a hardcoded member, summing it into the nodal acceleration.
- Time-varying or spatially varying loads can then be implemented as separate packages that write into `"bodyForce"` before `FEMLinearElasticity` runs — same multi-package composition pattern already used elsewhere.

The challenge is ordering: `FEMLinearElasticity` needs the body force field to be populated before its own `EvaluateDerivatives`. The integrator calls packages in order, so the body force package should be listed first in the package vector. Alternatively, `FEMLinearElasticity` could call a virtual `computeBodyForce()` method that subclasses override.

## Files Changed / Created

```
src/Materials/
  material.hh            — FEMMaterial<dim> abstract base
  linearElastic.cc       — IsotropicLinearElastic<dim>
  material.py            — Python binding (no trampoline)
  linearElastic.py       — Python binding
  Materials_PYB11.py     — module entry point
  CMakeLists.txt         — PYB11Generator_add_module(Materials)

src/Mesh/element.hh      — added computeStructuralStiffnessMatrix, computeStrain;
                           patched computeLumpedMassMatrix, computeStiffnessMatrix
                           to use std::abs(detJ) throughout QuadElement

src/Physics/
  femLinearElasticity.cc  — new physics package
  femLinearElasticity.py  — Python binding
  Physics_PYB11.py        — added Materials includes + femLinearElasticity imports

CMakeLists.txt (root)    — add_subdirectory(src/Materials)

tests/
  fem_elasticity_test.py  — cantilever test with SiloDump output
  yggdrasil.py            — added "Materials" to build path list
```
