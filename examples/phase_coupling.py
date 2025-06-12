from yggdrasil import *
from numpy import pi,sqrt
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

from matplotlib.colors import LinearSegmentedColormap
colors = [(1,0,0), (0, 0, 0), (0,0,1)]  # Red -> Black -> Blue
cmap = LinearSegmentedColormap.from_list('rbbl', colors, N=256)

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                       nx = 50,
                                       ny = 50,
                                       couplingConstant = 0.1,
                                       searchRadius = 5,
                                       dtmin = 0.001,
                                       save_as=None)
    
    constants = MKS()
    numNodes = nx*ny
    myNodeList = NodeList(numNodes)
    mesh = Grid2d(nx,ny,1.0/nx,1.0/ny)

    searchRadius = searchRadius/sqrt(nx*nx + ny*ny)

    phaseCoupling = PhaseCoupling2d(myNodeList, constants, 
                                couplingConstant = couplingConstant,
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

    integrator  = RungeKutta2Integrator2d(packages=packages, dtmin=dtmin, verbose=False)

    controller = Controller(integrator=integrator, periodicWork=[], statStep=10)

    getter = GetStrength(strength,mesh)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=getter,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="phi")
        AnimateGrid2d(bounds,update_method,extremis=[-1,1],cmap=cmap,frames=400,save_as=save_as)
    else:
        controller.Step(5000)
