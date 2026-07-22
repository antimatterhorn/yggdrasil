# Copyright (C) 2026  Cody Raskin

from LinearAlgebra import *
from numpy import pi, cos, sin


class ConstantDThetaPolyGenerator:
    """Like ConstantNThetaPolyGenerator (also exposes `.cells`, not just
    `.positions`, unlike this module's other generators), but ntheta grows
    outward instead of staying constant on every ring -- doubling whenever
    the constant-arc-length target (the same one ConstantDThetaDisk2d uses)
    has grown to at least 2x the current ring's ntheta -- so cells don't get
    wide angularly at large radius the way ConstantNThetaPolyGenerator's do.
    Consecutive rings connect with plain quads where ntheta is unchanged,
    and with the standard 1-inner-segment-to-2-outer-segments triangle
    pattern (3 triangles, sharing the new outer midpoint vertex) wherever
    ntheta doubles. The single center point is still a triangle fan to the
    innermost ring, same as ConstantNThetaPolyGenerator.
    """
    def __init__(self, nrings, ntheta0, thetaMin=0.0, thetaMax=2 * pi):
        self.positions = []
        self.cells = []

        dtheta = thetaMax - thetaMin
        fullCircle = abs(dtheta - 2 * pi) < 1e-12
        # Target angular spacing, calibrated so ntheta0 is reached at the
        # innermost ring (radius 1/nrings), not at the outer radius 1 --
        # d = dtheta/ntheta0 would put ntheta0 at r=1 instead, so the ideal
        # ntheta could never exceed ntheta0 anywhere inside the disk and
        # doubling would never trigger.
        d = dtheta / (nrings * ntheta0)

        self.positions.append([0.0, 0.0])  # center, node 0

        ringStart = []   # ringStart[k] = node index of theta=thetaMin on ring k+1
        ringNtheta = []  # ringNtheta[k] = ntheta on ring k+1

        ntheta = ntheta0
        nextIndex = 1
        for i in range(1, nrings + 1):
            r = i / float(nrings)
            idealNtheta = max(1, round(dtheta * r / d))
            if idealNtheta >= 2 * ntheta:
                ntheta *= 2

            pointsPerRing = ntheta if fullCircle else ntheta + 1
            ringStart.append(nextIndex)
            ringNtheta.append(ntheta)
            for j in range(pointsPerRing):
                theta = thetaMin + dtheta * j / ntheta
                self.positions.append([r * cos(theta), r * sin(theta)])
            nextIndex += pointsPerRing

        def ringNode(i, j):
            # i is the 1-based ring index; j wraps for a full circle.
            n = ringNtheta[i - 1]
            pointsPerRing = n if fullCircle else n + 1
            return ringStart[i - 1] + (j % pointsPerRing)

        # Triangle fan from the center to the innermost ring.
        for j in range(ringNtheta[0]):
            self.cells.append([0, ringNode(1, j), ringNode(1, j + 1)])

        # Ring i to ring i+1: plain quads, or a doubling transition.
        for i in range(1, nrings):
            nInner = ringNtheta[i - 1]
            nOuter = ringNtheta[i]
            if nOuter == nInner:
                for j in range(nInner):
                    self.cells.append([ringNode(i, j), ringNode(i, j + 1),
                                        ringNode(i + 1, j + 1), ringNode(i + 1, j)])
            else:
                for j in range(nInner):
                    a0, a1 = ringNode(i, j), ringNode(i, j + 1)
                    b0 = ringNode(i + 1, 2 * j)
                    b1 = ringNode(i + 1, 2 * j + 1)
                    b2 = ringNode(i + 1, 2 * j + 2)
                    self.cells.append([a0, b0, b1])
                    self.cells.append([a0, b1, a1])
                    self.cells.append([a1, b1, b2])

        self.npoints = len(self.positions)
        self.ncells = len(self.cells)
