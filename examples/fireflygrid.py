from yggdrasil import *
from numpy import pi
from random import random
from Physics import PhaseCoupling2d
from Mesh import Grid2d
from Animation import *

class GetStrength:
    def __init__(self,sField,mesh):
        self.field = sField
        self.mesh = mesh
    def __call__(self,i,j,fieldName):
        id = self.mesh.index(i,j,0)
        return self.field[id]

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                       nx = 50,
                                       ny = 50,
                                       couplingConstant = 0.1,
                                       lightFraction = 0,
                                       searchRadius = 0,
                                       dtmin = 0.01)
    
    constants = MKS()
    numNodes = nx*ny
    myNodeList = NodeList(numNodes)
    mesh = Grid2d(nx,ny,1,1)
    phaseCoupling = PhaseCoupling2d(myNodeList, constants, 
                                couplingConstant = couplingConstant, lightFraction = lightFraction,
                                searchRadius = searchRadius)
    packages = [phaseCoupling]

    mesh.assignPositions(myNodeList)

    phase = myNodeList.getFieldDouble("kphase")
    pos = myNodeList.getFieldVector2d("position")
    strength = myNodeList.getFieldDouble("kstrength")
    omega = myNodeList.getFieldDouble("komega")
    for i in range(numNodes):
        phase.setValue(i,random()*2*pi)
        omega.setValue(i,random()*0.1+1.0)

    integrator  = Integrator2d(packages=packages, dtmin=dtmin, verbose=False)

    controller = Controller(integrator=integrator, periodicWork=[], statStep=10)

    getter = GetStrength(strength,mesh)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=getter,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="phi")
        AnimateGrid2d(bounds,update_method,extremis=[-1,1],cmap="plasma",frames=1000)
    else:
        controller.Step(5000)
