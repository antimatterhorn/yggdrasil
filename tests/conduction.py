from yggdrasil import *
import matplotlib.pyplot as plt
from Animation import *


if __name__ == "__main__":
    commandLine = CommandLineArguments(animate = True,
                                        siloDump = False,
                                        cycles = 3000,
                                        nx = 100,
                                        ny = 100,
                                        dx = 1,
                                        dy = 1,
                                        dtmin = 0.001,
                                        intVerbose = False)

    myGrid = Grid2d(nx,ny,dx,dy)
    print("grid size:",myGrid.size())
    
    myNodeList = NodeList(nx*ny)
    print("numNodes =",myNodeList.numNodes)
    print("field names =",myNodeList.fieldNames)

    constants = MKS()
    eos = IdealGasEOS(1.4,constants)
    print(eos,"gamma =",eos.gamma)
    opac = ConstantOpacity(0.15,constants)

    cond = ThermalConduction2d(myNodeList,constants,eos,opac,myGrid)

    box = ReflectingGridBoundaries2d(grid=myGrid)
    cond.addBoundary(box)

    integrator = RungeKutta4Integrator2d([cond],dtmin=dtmin,verbose=intVerbose)

    density = myNodeList.getFieldDouble("density")
    energy  = myNodeList.getFieldDouble("specificInternalEnergy")
    position = myNodeList.getFieldVector2d("position")

    p0 = 2.5
    gamma = eos.gamma

    for j in range(ny):
        for i in range(nx):
            idx = myGrid.index(i,j,0)
            energy.setValue(idx, 20.0)
            density.setValue(idx, 1.0)
            pos = position[idx]
            x = pos.x - nx//2
            y = pos.y - ny//2
            r = np.sqrt(x*x + y*y)
            if r < nx//6:
                energy.setValue(idx, 50.0)



    periodicWork = []

    if siloDump:
        meshWriter = SiloDump(baseName="cond",
                                nodeList=myNodeList,
                                fieldNames=["density","specificInternalEnergy","pressure"],
                                dumpCycle=50)
        periodicWork += [meshWriter]

    controller = Controller(integrator=integrator,periodicWork=periodicWork,statStep=50)

    if(animate):
        title = MakeTitle(controller,"time","time")

        bounds = (nx,ny)
        update_method = AnimationUpdateMethod2d(call=cond.getCell2d,
                                                stepper=controller.Step,
                                                title=title,
                                                fieldName="temperature")
        AnimateGrid2d(bounds,update_method,extremis=[0,5],frames=cycles,cmap="plasma")
    else:
        controller.Step(cycles)
