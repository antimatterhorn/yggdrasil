from yggdrasil import *
from Animation import AnimateScatter
from Physics import ConstantForce2d,Kinetics2d

if __name__ == "__main__":
    constants = MKS()
    numNodes = 1
    myNodeList = NodeList(numNodes)
    # Node_created
    gravVec = Vector2d(0, -9.8)

    constantForce   = ConstantForce2d(myNodeList, constants, gravVec)
    kinetics        = Kinetics2d(myNodeList,constants)
    packages        = [constantForce,kinetics]
    # Packages_created
    integrator  = RungeKutta4Integrator2d(packages=packages, dtmin=0.01, verbose=False)
    print(integrator)
    # Integrator_created
    rad     = myNodeList.radius
    mass    = myNodeList.mass
    vel     = myNodeList.velocity
    pos     = myNodeList.position
    
    mass[0] = 1
    rad[0] = 0.01
    vel[0] = Vector2d(5,5)
    pos[0] = Vector2d(0,0)
    # Initial_conditions_set
    controller = Controller(integrator=integrator, periodicWork=[], statStep=1)

    bounds = (-1, 6, -1, 1.4)
    AnimateScatter(bounds, stepper=controller, positions=pos, frames=10, interval=50)