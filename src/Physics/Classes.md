
```mermaid
classDiagram
    Kinematics <|-- PointSourceGravity
    Kinematics <|-- ConstantGravity
    Kinematics <|-- NBodyGravity
    Physics <|-- PhaseCoupling
    Physics <|-- WaveEquation
    Physics <|-- Hydro
    Physics <|-- Kinematics
    Kinematics <|-- Kinetics
    Physics : +NodeList* nodeList
    Physics : +PhysicalConstants& constants
    Physics : VerifyFields(NodeList* nodeList)
    Physics : ZeroTimeInitialize()
    Physics : PrestepInitialize()
    Physics : EvaluateDerivatives(State* initialState, State<dim>& deriv, double time, double dt)
    Physics : FinalizeStep(State* finalState)
    Physics : FinalChecks()
    Physics : EnrollFields[typename T](string[] fields)
    Physics : EnrollStateFields[typename T](string[] fields)
    class PointSourceGravity{
        +Vector pointSourceLocation
        +Vector pointSourceVelocity
        +double pointSourceMass
    }
    class ConstantGravity{
        +Vector gravity
    }
    class NBodyGravity{
        +double plummerLength
    }
    class Hydro{
        +EquationOfState* eos
    }
    class WaveEquation{
        +Grid* grid
        +double soundSpeed
    }
    WaveEquation o -- _WaveEquation
    class _WaveEquation{
        +Grid* grid
        +string depthMap
    }
    class RockPaperScissors{
        +Grid* grid
        +double A
        +double D
    }
    Hydro <| -- GridHydro
    class GridHydro{
        +Grid* grid
    }
    GridHydro <| -- GridHydroHLLE
    class GridHydroHLLE{
    }
    GridHydro <| -- GridHydroHLLC
    class GridHydroHLLC{
    }
    Hydro <| -- EulerHydro
    class EulerHydro{
        +Grid* grid
    }
    class PhaseCoupling{
        +double couplingConstant
        [+double searchRadius]*
    }
```