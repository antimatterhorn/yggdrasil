```mermaid
classDiagram
Integrator <|-- RungeKutta2Integrator
Integrator <|-- RungeKutta4Integrator
Integrator <|-- CrankNicolsonIntegrator
Integrator : +Physics* physics
Integrator : +double dtmin
Integrator : Step()
Integrator : double Time()
Integrator : int Cycle()
Integrator : double Dt()
```