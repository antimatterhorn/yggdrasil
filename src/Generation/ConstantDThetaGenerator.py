# Copyright (C) 2026  Cody Raskin

from LinearAlgebra import *
from numpy import pi,sqrt,cos,sin

    
class ConstantDThetaDisk2d:
    def __init__(self, approx_npoints, thetaMin=0.0, thetaMax=2 * pi):
        self.positions = []
        d = sqrt(pi / approx_npoints)  # target spacing
        dr = d
        nrings = int(1.0 / dr + 0.5)
        dtheta = thetaMax - thetaMin
        fullCircle = abs(dtheta - 2 * pi) < 1e-12

        self.positions.append([0.0, 0.0])  # center

        for i in range(1, nrings + 1):
            r = i * dr
            if i == nrings:
                r = 1.0  # force outermost ring on unit circle

            ntheta = max(1, round(dtheta * r / d))

            # A full circle's last point coincides with its first, so it's
            # excluded; a partial sector needs both bounding edges included
            # so ring points land exactly on the sector's straight edges.
            jmax = ntheta if not fullCircle else ntheta - 1
            for j in range(jmax + 1):
                theta = thetaMin + dtheta * j / ntheta
                x = r * cos(theta)
                y = r * sin(theta)
                self.positions.append([x, y])

        # Optional: store actual number of points
        self.npoints = len(self.positions)
