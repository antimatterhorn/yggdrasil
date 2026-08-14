from yggdrasil import *
from Animation import *
from math import sin,cos
from Mesh import Grid2d
from Physics import GridHydroKT2d
from EOS import IdealGasEOS
from Boundaries import OutflowGridBoundary2d, DirichletGridBoundary2d

class oscillate:
    def __init__(self,nodeList,grid,width,height,workCycle=1):
        self.nodeList = nodeList
        self.cycle = workCycle
        self.grid = grid
        self.width = width
        self.height = height
        self.pressure = nodeList.pressure
    def __call__(self,cycle,time,dt):
        a = 5*(cos(time))
        i = int(self.width/6)
        j = int(self.height/4)
        idx = self.grid.index(i,j,0)
        self.pressure.setValue(idx,a)
        a = 5*(cos(time*1.1))
        i = int(self.width/6)
        j = int(self.height*3/4)
        idx = self.grid.index(i,j,0)
        self.pressure.setValue(idx,a)

from Utilities import SiloDump
      

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                       cycles = 20000,
                                       nx = 100,
                                       ny = 100)
    

    constants = MKS()

    grid = Grid2d(nx,ny,1,1)

    myNodeList = NodeList(grid.size())
    print("grid %dx%d"%(nx,ny))
    print(grid)

    eos = IdealGasEOS(1.4,constants)
    hydro = GridHydroKT2d(myNodeList,constants,eos,grid) 

    print(hydro)

    packages = [hydro]

    pm = OutflowGridBoundary2d(grid=grid)
    box = DirichletGridBoundary2d(grid=grid)

    nbox = 10
    dy = ny/nbox
    h  = int(dy/2.5)
    for i in range(nbox):
        x = int(nx/3)
        y = dy*(i+1)
        box.addSphere(Vector2d(x,y),2)

    hydro.addBoundary(pm)
    hydro.addBoundary(box)

    integrator = RungeKutta4Integrator2d(packages=packages,
                              dtmin=1e-4,verbose=False)
    print(integrator)

    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    # INITIAL CONDITIONS

    density = myNodeList.getFieldDouble("density")
    energy  = myNodeList.getFieldDouble("specificInternalEnergy")
    velocity = myNodeList.getFieldVector2d("velocity")
    position = myNodeList.getFieldVector2d("position")

    p0 = 2.5
    p_src = 20.0   # elevated pressure inside source circles
    r_src = 6      # source circle radius in cells
    gamma = eos.gamma

    sources = [(nx // 6, ny // 4), (nx // 6, 3 * ny // 4)]

    for j in range(ny):
        for i in range(nx):
            idx = grid.index(i, j, 0)
            rho = 1.0
            p = p0

            for (cx, cy) in sources:
                if (i - cx)**2 + (j - cy)**2 <= r_src**2:
                    p = p_src
                    break

            density.setValue(idx, rho)
            energy.setValue(idx, p / ((gamma - 1.0) * rho))
            velocity.setValue(idx, Vector2d(0.0, 0.0))

    #osc = oscillate(nodeList=myNodeList,grid=grid,width=nx,height=ny,workCycle=1)
    periodicWork = []
    if (not animate):
        vtk = SiloDump("testMesh",myNodeList,fieldNames=["density","pressure","velocity"],dumpCycle=50,grid=grid)
        periodicWork.append(vtk)

    controller = Controller(integrator=integrator,
                            statStep=100,
                            periodicWork=periodicWork)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=hydro.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="pressure") # change to 'phi' to view full wave
        AnimateGrid2d(bounds,update_method,extremis=[0,25],frames=cycles,cmap='plasma')
    else:
        controller.Step(cycles)
