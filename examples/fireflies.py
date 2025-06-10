from yggdrasil import *
from numpy import pi
from random import random
from Physics import PhaseCoupling2d, Kinetics2d
from RandomNodeGenerator import RandomNodeGenerator2d
from Animation import AnimateScatter

if __name__ == "__main__":
    commandLine = CommandLineArguments(numNodes = 100,
                                       couplingConstant = 0.1,
                                       lightFraction = 0.2)
    
    constants = MKS()
    myNodeList = NodeList(numNodes)

    phaseCoupling = PhaseCoupling2d(myNodeList, constants, 
                                    couplingConstant = couplingConstant, lightFraction = lightFraction)
    kinetics = Kinetics2d(myNodeList,constants)
    packages = [phaseCoupling,kinetics]

    myNodeList.updatePositions2d(RandomNodeGenerator2d(numNodes).positions)

    rad = myNodeList.getFieldDouble("radius")
    mass = myNodeList.getFieldDouble("mass")
    phase = myNodeList.getFieldDouble("kphase")
    pos = myNodeList.getFieldVector2d("position")
    strength = myNodeList.getFieldDouble("kstrength")
    omega = myNodeList.getFieldDouble("komega")
    for i in range(numNodes):
        rad.setValue(i, 1e-5)
        mass.setValue(i,1e-5)
        phase.setValue(i,random()*2*pi)
        omega.setValue(i,random()*0.1+1.0)

    integrator  = Integrator2d(packages=packages, dtmin=0.1, verbose=False)

    controller = Controller(integrator=integrator, periodicWork=[], statStep=1000)

    bounds = (-1, 1, -1, 1)
    AnimateScatter(bounds, stepper=controller, positions=pos, 
                   get_color_field=lambda i: strength[i],
                   cmap='plasma',
                   color_limits=(0.0, 1.0),
                   background="black",  # or "#202020" or (0.1, 0.1, 0.1))
                    )
