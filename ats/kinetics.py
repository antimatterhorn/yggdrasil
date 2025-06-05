from yggdrasil import *
from Animation import *

def run():
    animate = False
    cycles = 35
    numNodes = 2

    constants = MKS()
    myNodeList = NodeList(numNodes)

    gravVec = Vector2d(0, 0)

    constantGravity = ConstantGravity2d(myNodeList, constants, gravVec)
    kinetics = Kinetics2d(myNodeList,constants)
    packages = [constantGravity,kinetics]

    integrator = RungeKutta2Integrator2d(packages=packages, dtmin=0.01,verbose=False)


    rad = myNodeList.getFieldDouble("radius")
    mass = myNodeList.getFieldDouble("mass")
    for i in range(numNodes):
        rad.setValue(i,2)
        mass.setValue(i,0.2)

    pos = myNodeList.getFieldVector2d("position")
    vel = myNodeList.getFieldVector2d("velocity")
    for i in range(numNodes):
        pos.setValue(i, Vector2d(-10+i/numNodes*30, -0.5+i/numNodes*1))
        vel.setValue(i, Vector2d(3*(-1)**i,0))

    controller = Controller(integrator=integrator, periodicWork=[], statStep=100)

    bounds = (-10, 10, -10, 10)
    if animate:
        AnimateScatter(bounds, stepper=controller, positions=pos, frames=10, interval=50)
    else:
        out = []
        controller.Step(cycles)
        out.append(integrator.time)
        for i in range(numNodes):
            out.append(pos[i].x)
        return out

if __name__ == "__main__":
    print(run())
