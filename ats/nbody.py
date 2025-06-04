from yggdrasil import *
from Animation import *


class dumpState:
    def __init__(self,nodeList,workCycle=1,G=1):
        self.nodeList = nodeList
        self.cycle = workCycle
        self.dump = []
        self.energy = []
        self.G = G
    def __call__(self,cycle,time,dt):
        for i in range(self.nodeList.numNodes):
            self.dump.append((self.nodeList.getFieldVector2d("position")[i].x,self.nodeList.getFieldVector2d("position")[i].y))


def run():
    numNodes = 3
    cycles = 330
    animate = False
    myNodeList = NodeList(numNodes)

    constants = PhysicalConstants(6.378e+6,     # earth radius in m
                                  5.972e+24,    # earth mass in kg 
                                  1.0,          # s
                                  1.0, 
                                  1.0) 
    loc = Vector2d(0, 0)

    nBodyGrav = NBodyGravity2d(nodeList=myNodeList,
                               constants=constants,
                               plummerLength=0.01)
    packages = [nBodyGrav]

    positions   = myNodeList.getFieldVector2d("position")
    mass        = myNodeList.getFieldDouble("mass")
    for i in range(numNodes):
        mass.setValue(i,1.0)

    positions.setValue(0,Vector2d(-1,1))
    positions.setValue(1,Vector2d(1,1))
    positions.setValue(2,Vector2d(0,-1))

    integrator = RungeKutta4Integrator2d(packages=packages,
                                         dtmin=0.5e1,verbose=False)
  

    dump = dumpState(myNodeList,workCycle=1000,G=constants.G)
    periodicWork = [dump]

    controller = Controller(integrator=integrator,periodicWork=[],statStep=499)

    bounds = (-1,1,-1,1)
    if animate:
        AnimateScatter(bounds, stepper=controller, positions=positions, frames=100, interval=50)
    else:
        out = []
        controller.Step(cycles)
        out.append(integrator.Time())
        for i in range(numNodes):
            out.append(positions[i].y)
        return out

if __name__ == "__main__":
    print(run())
