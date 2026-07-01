# src

C++ core of the Yggdrasil physics simulation framework, exposed to Python via pybind11/PYB11Generator. Organized into a foundation layer (Type, Math, DataBase, State) with no circular dependencies, a simulation layer built on top of it (Mesh, EOS, Materials, Trees, Boundaries, Physics, Integrators), and Python-only runtime support (Generation, IO, Utilities, Calculators).

## Contents

| Directory | Description |
|---|---|
| [Type/](Type/) | Foundational types everything depends on: `PhysicalConstants` (unit system, always passed by reference) and `Name` (field identifier) |
| [Math/](Math/) | `Lin::Vector<dim>` and `Lin::Tensor<dim>` linear algebra primitives, plus slope limiters |
| [DataBase/](DataBase/) | `Field<T>` (named typed array) and `NodeList` (field container) — the core data structures everything else builds on |
| [State/](State/) | `State<dim>`, an integration-step snapshot of enrolled fields passed between `Physics` and `Integrator` |
| [Mesh/](Mesh/) | Geometry/topology: `Grid` (structured Cartesian), `FEMesh` (unstructured FE), `VoronoiMesh` (tessellation), and `Element` types |
| [EOS/](EOS/) | Equation-of-state and opacity models (ideal gas, isothermal, etc.) providing thermodynamic closure to hydro solvers |
| [Materials/](Materials/) | FEM structural constitutive models (`FEMMaterial`, `LinearElastic`) producing the material matrix D for stiffness assembly |
| [Trees/](Trees/) | `KDTree` for neighbor search and `SpatialTree` (Barnes-Hut) for hierarchical N-body gravity |
| [Boundaries/](Boundaries/) | `Boundary`, `GridBoundary` variants (reflecting, periodic, outflow, Dirichlet), `Collider`, and `Constraint` |
| [Physics/](Physics/) | All physics packages derived from `Physics<dim>`: hydro, gravity, FEM, wave equation, kinetics, magnetic field, CNC pathing, and more |
| [Integrators/](Integrators/) | `Integrator<dim>` base and concrete time integrators: RK2, RK4, Crank-Nicolson |
| [Generation/](Generation/) | Pure-Python node/particle position generators (lattice, HCP, Poisson, Voronoi/CVT, tanh-stretched, etc.) |
| [IO/](IO/) | `SiloMeshWriter`, `VtkMeshWriter`, and importers for depth maps and OBJ geometry |
| [Utilities/](Utilities/) | Python-side runtime support: `Controller` (time loop driver), animation, CLI args, progress bar |
| [Calculators/](Calculators/) | Standalone tools not integrated into the Physics/Integrator framework: `Cosmology`, `Mandelbrot`, `StringArt`, `TimeDilation` |
