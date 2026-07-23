# Copyright (C) 2026  Cody Raskin

from LinearAlgebra import *
import numpy as np
from scipy.spatial import Delaunay
from ConstantDThetaGenerator import ConstantDThetaDisk2d
from PointsToPoly import PointsToPoly2d

class ConstantDThetaPolyDisk2d:
    def __init__(self, npoints, thetaMin=0.0, thetaMax=2 * np.pi, angleTol=25.0):
        self.positions = ConstantDThetaDisk2d(npoints, thetaMin, thetaMax).positions

        self.poly = PointsToPoly2d(self.positions, angleTol)
        self.cells = self.poly.cells

        self.npoints = len(self.positions)
        self.ncells = len(self.cells)
