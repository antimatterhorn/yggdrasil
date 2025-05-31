from LinearAlgebra import *
from numpy import pi,sqrt,cos,sin

    
class ConstantDThetaDisk2d:
    def __init__(self, approx_npoints):
        self.positions = []
        d = sqrt(pi / approx_npoints)  # target spacing
        dr = d
        nrings = int(1.0 / dr + 0.5)

        self.positions.append([0.0, 0.0])  # center

        for i in range(1, nrings + 1):
            r = i * dr
            if i == nrings:
                r = 1.0  # force outermost ring on unit circle

            ntheta = max(1, round(2 * pi * r / d))

            for j in range(ntheta):
                theta = 2 * pi * j / ntheta
                x = r * cos(theta)
                y = r * sin(theta)
                self.positions.append([x, y])

        # Optional: store actual number of points
        self.npoints = len(self.positions)
