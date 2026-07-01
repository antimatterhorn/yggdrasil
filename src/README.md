# src

C++ core of the Yggdrasil physics simulation framework, exposed to Python via pybind11/PYB11Generator. Organized into a foundation layer (Type, Math, DataBase, State) with no circular dependencies, a simulation layer built on top of it (Mesh, EOS, Materials, Trees, Boundaries, Physics, Integrators), and Python-only runtime support (Generation, IO, Utilities, Calculators). See the top-level [CLAUDE.md](../CLAUDE.md) for the full architecture diagram; each subdirectory below has its own `CLAUDE.md` with class-level documentation.

## Contents

| Directory | Description |
|---|---|
| [Type/](Type/CLAUDE.md) | Foundational types everything depends on: `PhysicalConstants` (unit system, always passed by reference) and `Name` (field identifier) |
| [Math/](Math/CLAUDE.md) | `Lin::Vector<dim>` and `Lin::Tensor<dim>` linear algebra primitives, plus slope limiters |
| [DataBase/](DataBase/CLAUDE.md) | `Field<T>` (named typed array) and `NodeList` (field container) — the core data structures everything else builds on |
| [State/](State/CLAUDE.md) | `State<dim>`, an integration-step snapshot of enrolled fields passed between `Physics` and `Integrator` |
| [Mesh/](Mesh/CLAUDE.md) | Geometry/topology: `Grid` (structured Cartesian), `FEMesh` (unstructured FE), `VoronoiMesh` (tessellation), and `Element` types |
| [EOS/](EOS/CLAUDE.md) | Equation-of-state and opacity models (ideal gas, isothermal, etc.) providing thermodynamic closure to hydro solvers |
| [Materials/](Materials/CLAUDE.md) | FEM structural constitutive models (`FEMMaterial`, `LinearElastic`) producing the material matrix D for stiffness assembly |
| [Trees/](Trees/CLAUDE.md) | `KDTree` for neighbor search and `SpatialTree` (Barnes-Hut) for hierarchical N-body gravity |
| [Boundaries/](Boundaries/CLAUDE.md) | `Boundary`, `GridBoundary` variants (reflecting, periodic, outflow, Dirichlet), `Collider`, and `Constraint` |
| [Physics/](Physics/CLAUDE.md) | All physics packages derived from `Physics<dim>`: hydro, gravity, FEM, wave equation, kinetics, magnetic field, CNC pathing, and more |
| [Integrators/](Integrators/CLAUDE.md) | `Integrator<dim>` base and concrete time integrators: RK2, RK4, Crank-Nicolson |
| [Generation/](Generation/CLAUDE.md) | Pure-Python node/particle position generators (lattice, HCP, Poisson, Voronoi/CVT, tanh-stretched, etc.) |
| [IO/](IO/CLAUDE.md) | `SiloMeshWriter`, `VtkMeshWriter`, and importers for depth maps and OBJ geometry |
| [Utilities/](Utilities/CLAUDE.md) | Python-side runtime support: `Controller` (time loop driver), animation, CLI args, progress bar |
| [Calculators/](Calculators/CLAUDE.md) | Standalone tools not integrated into the Physics/Integrator framework: `Cosmology`, `Mandelbrot`, `StringArt`, `TimeDilation` |
