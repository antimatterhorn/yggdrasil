from yggdrasil import *
from Animation import *
from math import sin,cos
from Mesh import Grid2d
from Physics import WaveEquation2d
from Boundaries import OutflowGridBoundary2d, DirichletGridBoundary2d

class oscillate:
    def __init__(self,nodeList,grid,width,height,workCycle=1):
        self.nodeList = nodeList
        self.cycle = workCycle
        self.grid = grid
        self.width = width
        self.height = height
        self.phi = myNodeList.getFieldDouble("phi")
    def __call__(self,cycle,time,dt):
        a = 5*(cos(time))
        i = 2
        for j in range(self.height):
            idx = self.grid.index(i,j,0)
            self.phi.setValue(idx,a)

from Utilities import SiloDump
      

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                       save = None, # doubleSlit.gif
                                       cycles = 800,
                                       nx = 100,
                                       ny = 100)
    

    constants = MKS()

    myNodeList = NodeList(nx*ny)
    
    grid = Grid2d(nx,ny,1,1)
    print("grid %dx%d"%(nx,ny))
    print(grid)

    waveEqn = WaveEquation2d(nodeList=myNodeList,
                             constants=constants,
                             grid=grid,C=1.0)

    print(waveEqn)

    packages = [waveEqn]

    pm = OutflowGridBoundary2d(grid=grid)
    box = DirichletGridBoundary2d(grid=grid)

    x = nx//3
    dx = 2
    box.addBox(Vector2d(x-dx,0),Vector2d(x,ny))
    s1 = 0.55*ny
    s2 = 0.45*ny

    box.removeBox(Vector2d(x-dx,s2-1),Vector2d(x,s2+1))
    box.removeBox(Vector2d(x-dx,s1-1),Vector2d(x,s1+1))

    waveEqn.addBoundary(pm)
    waveEqn.addBoundary(box)

    integrator = RungeKutta4Integrator2d(packages=packages,
                              dtmin=0.05,verbose=False)
    print(integrator)

    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    osc = oscillate(nodeList=myNodeList,grid=grid,width=nx,height=ny,workCycle=1)
    periodicWork = [osc]
    if (not animate):
        vtk = SiloDump("testMesh",myNodeList,fieldNames=["phi","xi"],dumpCycle=50)
        periodicWork.append(vtk)

    controller = Controller(integrator=integrator,
                            statStep=100,
                            periodicWork=periodicWork)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=waveEqn.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="waveEnergyDensity") # change to 'phi' to view full wave
        AnimateGrid2d(bounds,update_method,extremis=[0,0.03],frames=cycles,cmap='plasma',save_as=save)
    else:
        controller.Step(cycles)
