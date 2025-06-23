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
        self.phi = self.nodeList.getFieldDouble("phi")
    def __call__(self,cycle,time,dt):
        a = 5*(cos(0.5*time))*np.exp(-time)
        i = int(self.width/2)
        j = int(self.height/2)
        idx = self.grid.index(i,j,0)
        self.phi.setValue(idx,a)

class SpkOscillator:
    def __init__(self,speaker,nodeList,grid,width,height,locations,workCycle=1):
        self.speaker = speaker
        self.nodeList = nodeList
        self.cycle = workCycle
        self.grid = grid
        self.width = width
        self.height = height
        self.phi = self.nodeList.getFieldDouble("phi")
        self.locations = locations
    def __call__(self,cycle,time,dt):
        m = self.speaker.get_output(cycle*2.5e-3)
        #print(cycle*2.5e-3,a)
        for loc in self.locations:
            i = int(loc[0])
            j = int(loc[1])
            idx = self.grid.index(i,j,0)
            self.phi.setValue(idx,m)

from Utilities import DampedHarmonicOscillator,SiloDump 

if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                       dump = False,
                                       cycles = 120000,
                                       nx = 150,
                                       ny = 80,
                                       cs = 1.0,
                                       frequency = 10.0,
                                       amplitude  = 8.0,
                                       damping = -5.0,
                                       ns = 16)
    
    constants = MKS()

    myNodeList = NodeList(nx*ny)
    
    grid = Grid2d(nx,ny,1,1)
    print("grid %dx%d"%(nx,ny))
    print(grid)

    waveEqn = WaveEquation2d(nodeList=myNodeList,
                             constants=constants,
                             grid=grid,C=cs)

    print(waveEqn)

    packages = [waveEqn]

    box = DirichletGridBoundary2d(grid=grid)

    box.addBox(Vector2d(12,32),Vector2d(16,48))
    box.addBox(Vector2d(12,20),Vector2d(16,28))
    box.addBox(Vector2d(16,48),Vector2d(20,52))
    box.addBox(Vector2d(24,48),Vector2d(28,52))
    box.addBox(Vector2d(16,28),Vector2d(20,32))
    box.addBox(Vector2d(20,40),Vector2d(24,48))
    box.addBox(Vector2d(16,16),Vector2d(28,20))
    box.addBox(Vector2d(28,20),Vector2d(32,48))
    box.addBox(Vector2d(44,20),Vector2d(48,48))
    box.addBox(Vector2d(60,20),Vector2d(64,48))
    box.addBox(Vector2d(76,28),Vector2d(80,56))
    box.addBox(Vector2d(60,28),Vector2d(88,32))
    box.addBox(Vector2d(92,28),Vector2d(144,32))
    box.addBox(Vector2d(76,48),Vector2d(92,52))
    box.addBox(Vector2d(96,48),Vector2d(128,52))
    box.addBox(Vector2d(32,28),Vector2d(36,32))
    box.addBox(Vector2d(48,28),Vector2d(52,32))
    box.addBox(Vector2d(32,16),Vector2d(44,20))
    box.addBox(Vector2d(48,16),Vector2d(60,20))
    box.addBox(Vector2d(32,48),Vector2d(44,52))
    box.addBox(Vector2d(48,48),Vector2d(60,52))
    box.addBox(Vector2d(36,36),Vector2d(40,44))
    box.addBox(Vector2d(52,36),Vector2d(56,44))
    box.addBox(Vector2d(68,36),Vector2d(72,44))
    box.addBox(Vector2d(92,32),Vector2d(96,48))
    box.addBox(Vector2d(108,32),Vector2d(112,48))
    box.addBox(Vector2d(124,32),Vector2d(128,56))
    box.addBox(Vector2d(132,32),Vector2d(136,56))
    box.addBox(Vector2d(140,32),Vector2d(144,56))
    box.addBox(Vector2d(84,32),Vector2d(88,36))
    box.addBox(Vector2d(88,36),Vector2d(92,40))
    box.addBox(Vector2d(100,40),Vector2d(104,44))
    box.addBox(Vector2d(112,36),Vector2d(116,40))
    box.addBox(Vector2d(120,40),Vector2d(124,44))
    box.addBox(Vector2d(128,44),Vector2d(132,48))
    box.addBox(Vector2d(128,56),Vector2d(132,60))
    box.addBox(Vector2d(136,56),Vector2d(140,60))
    box.addBox(Vector2d(72,56),Vector2d(76,60))
    box.addBox(Vector2d(68,48),Vector2d(72,56))
    box.addBox(Vector2d(64,48),Vector2d(72,52))



    box.addDomain()

    waveEqn.addBoundary(box)

    integrator = RungeKutta4Integrator2d(packages=packages,dtmin=0.05,verbose=False)
    print(integrator)

    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    locations = []
    for x in [20,36,52,84,100,116,132]:
        locations.append([x,64])
    for x in [128,136]:
        locations.append([x,52])
    for x in [65,81]:
        locations.append([x,47])
    for x in [116,129]:
        locations.append([x,44])
    locations.append([100,40])
    for x in [20,36,52,72]:
        locations.append([x,36])
    for x in [68,84,100,116,132]:
        locations.append([x,24])
    for x in [20,36,52]:
        locations.append([x,12])

    locations.append([8,64])
    locations.append([8,32])
    locations.append([8,12])
    locations.append([40,48])
    locations.append([56,48])
    locations.append([120,36])



    #osc = oscillate(nodeList=myNodeList,grid=grid,width=nx,height=ny,workCycle=1)
    spk = DampedHarmonicOscillator(frequency=frequency,amplitude=amplitude,damping=damping)
    osc = SpkOscillator(nodeList=myNodeList,grid=grid,width=nx,height=ny,locations=locations,workCycle=1,speaker=spk)
    periodicWork = [osc]
    if(dump):
        silo = SiloDump("testMesh",myNodeList,fieldNames=["phi","xi"],dumpCycle=50)
        periodicWork.append(silo)

    controller = Controller(integrator=integrator,
                            statStep=500,
                            periodicWork=periodicWork)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=waveEqn.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="maxphi")
        AnimateGrid2d(bounds,update_method,extremis=[-5,5],frames=cycles,cmap="plasma")
    else:
        controller.Step(cycles)
