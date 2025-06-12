from yggdrasil import *
from numpy import pi
from random import random
from Physics import PhaseCoupling2d, Kinetics2d
from PoissonNodeGenerator import PoissonNodeGenerator2d
from Animation import AnimateScatter
from matplotlib.colors import LinearSegmentedColormap

colors = [(0, 0, 0), (0.7, 1, 0)]  # Black -> Green
cmap = LinearSegmentedColormap.from_list('green_black', colors, N=256)

class light:
    def __init__(self,field,frac):
        self.field = field
        self.frac = frac
    def __call__(self,i):
        return int(self.field[i] > (1.0-self.frac))

if __name__ == "__main__":
    commandLine = CommandLineArguments(numNodes = 200,
                                       couplingConstant = 0.2,
                                       lightFraction = 0.2,
                                       statStep = 1000)
    
    constants = MKS()




    dist = PoissonNodeGenerator2d(numNodes).positions
    numNodes = len(dist)
    myNodeList = NodeList(numNodes)
    phaseCoupling = PhaseCoupling2d(myNodeList, constants, 
                                couplingConstant = couplingConstant)
    kinetics = Kinetics2d(myNodeList,constants)
    packages = [phaseCoupling,kinetics]

    myNodeList.updatePositions2d(dist)

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

    integrator  = Integrator2d(packages=packages, dtmin=.001, verbose=False)

    controller = Controller(integrator=integrator, periodicWork=[], statStep=statStep)

    lightFrac = light(strength,lightFraction)

    bounds = (-1, 1, -1, 1)
    AnimateScatter(bounds, stepper=controller, positions=pos, 
                   get_color_field=lambda i: lightFrac(i),
                   cmap=cmap,
                   color_limits=(0.0, 1.0),
                   background="black",  # or "#202020" or (0.1, 0.1, 0.1))
                    )
