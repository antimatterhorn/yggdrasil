# Copyright (C) 2026  Cody Raskin

from LinearAlgebra import *
from numpy import pi, cos, sin


class ConstantNThetaPolyGenerator:
    """Node positions plus quad-dominant mesh connectivity for a polar disk
    or angular sector, unlike this module's other generators (which produce
    positions only). Rings share the same number of angular divisions
    (`ntheta`), so consecutive rings connect directly into quads; the single
    center point can't be a quad vertex, so the innermost ring connects to
    it as a triangle fan instead. `cells` is a list of node-index lists (4
    entries per quad, 3 per center-fan triangle) ready for
    `Mesh::ALEMesh.addCell`; `positions` is the usual flat [x, y] list this
    module's other generators expose, on a unit disk (or sector) for the
    caller to rescale.
    """
    def __init__(self, nrings, ntheta, thetaMin=0.0, thetaMax=2 * pi):
        self.positions = []
        self.cells = []

        dtheta = thetaMax - thetaMin
        fullCircle = abs(dtheta - 2 * pi) < 1e-12
        pointsPerRing = ntheta if fullCircle else ntheta + 1

        self.positions.append([0.0, 0.0])  # center, node 0

        def ringNode(i, j):
            # node index of angular position j (wrapped for a full circle)
            # on ring i (i=1..nrings); center is handled separately.
            return 1 + (i - 1) * pointsPerRing + (j % pointsPerRing)

        for i in range(1, nrings + 1):
            r = i / float(nrings)
            for j in range(pointsPerRing):
                theta = thetaMin + dtheta * j / ntheta
                self.positions.append([r * cos(theta), r * sin(theta)])

        # Triangle fan from the center to the innermost ring.
        for j in range(ntheta):
            self.cells.append([0, ringNode(1, j), ringNode(1, j + 1)])

        # Quads between consecutive rings.
        for i in range(1, nrings):
            for j in range(ntheta):
                self.cells.append([ringNode(i, j), ringNode(i, j + 1),
                                    ringNode(i + 1, j + 1), ringNode(i + 1, j)])

        self.npoints = len(self.positions)
        self.ncells = len(self.cells)
