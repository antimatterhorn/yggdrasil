from LinearAlgebra import *
from math import cos,sin,acos,asin,sqrt
from numpy import pi

class FibonacciSurface3d:
    # generates a surface of particles on a unit sphere using PS
    def __init__(self,nshell):
        self.positions = []
        fib = 1.6180397959
        ai  = 2.0*pi*fib
        for i in range(nshell):
            t = i/(float)(nshell)
            inc = acos(-1+2*t)
            az = ai * i
            x = sin(inc)*cos(az)
            y = sin(inc)*sin(az)
            z = cos(inc)
            

            self.positions.append([x,y,z])
        
        return
    
class FibonacciDisk2d:
    def __init__(self, npoints, radius=1.0, on_circle=False):
        self.positions = []
        golden_angle = pi * (3 - sqrt(5))  # ~2.39996 radians

        for i in range(npoints):
            theta = i * golden_angle
            if on_circle:
                r = radius  # All points lie on the unit circle
            else:
                r = radius * sqrt(i / npoints)  # Uniform disk density

            x = r * cos(theta)
            y = r * sin(theta)

            self.positions.append([x, y])