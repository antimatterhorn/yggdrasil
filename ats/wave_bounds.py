from yggdrasil import *
from Animation import *
from math import sin,cos

class oscillate:
    def __init__(self,nodeList,grid,width,height,workCycle=1):
        self.nodeList = nodeList
        self.cycle = workCycle
        self.grid = grid
        self.width = width
        self.height = height
        self.phi = self.nodeList.getFieldDouble("phi")
    def __call__(self,cycle,time,dt):
        a = 5*(cos(time))
        i = int(self.width/6)
        j = int(self.height/4)
        idx = self.grid.index(i,j,0)
        self.phi.setValue(idx,a)
        a = 5*(cos(time*1.1))
        i = int(self.width/6)
        j = int(self.height*3/4)
        idx = self.grid.index(i,j,0)
        self.phi.setValue(idx,a)

def run():    
    cycles = 200
    animate = False
    nx = 100
    ny = 100

    constants = MKS()

    myNodeList = NodeList(nx*ny)
    
    grid = Grid2d(nx,ny,1,1)

    waveEqn = WaveEquation2d(nodeList=myNodeList,
                             constants=constants,
                             grid=grid,C=1.0)


    packages = [waveEqn]

    pm = OutflowGridBoundary2d(grid=grid)
    box = DirichletGridBoundary2d(grid=grid)

    nbox = 10
    dy = ny/nbox
    h  = int(dy/2.5)
    for i in range(nbox):
        x = int(nx/3)
        y = dy*(i+1)
        box.addSphere(Vector2d(x,y),2)

    waveEqn.addBoundary(pm)
    waveEqn.addBoundary(box)

    integrator = RungeKutta4Integrator2d(packages=packages,
                              dtmin=0.05,verbose=False)

    osc = oscillate(nodeList=myNodeList,grid=grid,width=nx,height=ny,workCycle=1)
    periodicWork = [osc]

    controller = Controller(integrator=integrator,
                            statStep=1000,
                            periodicWork=periodicWork)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=waveEqn.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="maxphi") # change to 'phi' to view full wave
        AnimateGrid2d(bounds,update_method,extremis=[0,2],frames=cycles,cmap='plasma')
    else:
        out = []
        controller.Step(cycles)
        out.append(integrator.time)
        rgb_grid = np.zeros((ny, nx))
        max_values = []
        for j in range(ny):
            maxi = 0
            for i in range(nx):
                rgb_grid[j, i] = waveEqn.getCell2d(i % nx, j % ny, "maxphi")
                if i == nx // 2:
                    maxi = float(rgb_grid[j, i])
            max_values.append(maxi)
        return max_values


if __name__ == "__main__":
    print(run())