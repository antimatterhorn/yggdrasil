from yggdrasil import *
from Animation import *
from math import sin,cos
from Mesh import Grid2d
from Physics import WaveEquation2d,ComplexWaveEquation2d
from Boundaries import OutflowGridBoundary2d, DirichletGridBoundary2d

class oscillate:
    def __init__(self, nodeList, grid, width, height, workCycle=1):
        self.nodeList = nodeList
        self.cycle = workCycle
        self.grid = grid
        self.width = width
        self.height = height
        self.psi = self.nodeList.psi  # Updated field accessor

    def __call__(self, cycle, time, dt):
        a1 = 5 * cos(time)
        a2 = 5 * cos(1.1 * time)

        idx1 = self.grid.index(int(self.width / 6), int(self.height / 4), 0)
        idx2 = self.grid.index(int(self.width / 6), int(self.height * 3 / 4), 0)

        self.psi.setValue(idx1, complex(a1, 0.0))  # Inject real component
        self.psi.setValue(idx2, complex(0.0, a2))  # Inject imaginary component (for asymmetry)

from Utilities import SiloDump
      

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                       cycles = 20000,
                                       nx = 100,
                                       ny = 100)
    

    constants = MKS()

    myNodeList = NodeList(nx*ny)
    
    grid = Grid2d(nx,ny,1,1)
    print("grid %dx%d"%(nx,ny))
    print(grid)
    
    cwaveEqn = ComplexWaveEquation2d(nodeList=myNodeList,
                                    constants=constants,
                                    grid=grid,C=1.0)


    pm = OutflowGridBoundary2d(grid=grid)
    box = DirichletGridBoundary2d(grid=grid)

    nbox = 10
    dy = ny/nbox
    h  = int(dy/2.5)
    for i in range(nbox):
        x = int(nx/3)
        y = dy*(i+1)
        box.addSphere(Vector2d(x,y),2)

    cwaveEqn.addBoundary(pm)
    cwaveEqn.addBoundary(box)
    


    packages = [cwaveEqn]

    integrator = RungeKutta4Integrator2d(packages=packages,
                              dtmin=0.05,verbose=False)


    osc = oscillate(nodeList=myNodeList,grid=grid,width=nx,height=ny,workCycle=1)
    periodicWork = [osc]

    if (not animate):
        vtk = SiloDump("testMesh",myNodeList,fieldNames=["psi"],dumpCycle=50,grid=grid)
        periodicWork.append(vtk)

    controller = Controller(integrator=integrator,
                            statStep=100,
                            periodicWork=periodicWork)
    
    print(cwaveEqn)
    print(integrator)

    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=cwaveEqn.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="psi") # change to 'phi' to view full wave
        AnimateGrid2d(bounds,update_method,extremis=[-4,4],frames=cycles,cmap='plasma')
    else:
        controller.Step(cycles)
