from typing import Callable, Any
from yggdrasil import *
import numpy as np
from Utilities import MeshViz

def CircleQuad(n1,n2):
    from CircleQuadPolyGenerator import CircleQuadPolyGenerator2d
    return CircleQuadPolyGenerator2d(n1, n2, k=1.0)

def ConstantDTheta(n1,n2):
    from ConstantDThetaPolyGenerator import ConstantDThetaPolyDisk2d
    return ConstantDThetaPolyDisk2d(n1*n2)

def ConstantNTheta(n1,n2):
    from ConstantNThetaPolyGenerator import ConstantNThetaPolyGenerator
    return ConstantNThetaPolyGenerator(n1,n2)

Func = Callable[[int], Any]

getMesh: dict[str,Func] ={
    "CircleQuad": CircleQuad,
    "ConstantDTheta": ConstantDTheta,
    "ConstantNTheta": ConstantNTheta
}

if __name__ == "__main__":
    commandLine = CommandLineArguments(method="CircleQuad", n1 = 20, n2 = 20)

    gen = getMesh[method](n1, n2)

    mesh_viz = MeshViz(gen)
    mesh_viz.plot(method)