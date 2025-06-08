from yggdrasil import *
from Mesh import VoronoiMesh2d
from progressBar import ProgressBar
import math

class CVTNodeGenerator2d:
    def __init__(self, numPoints, n_iterations=10):
        self.positions = []

        # Step 1: Generate initial seed using Poisson Disk sampling
        from ConstantDThetaGenerator import ConstantDThetaDisk2d
        initial = ConstantDThetaDisk2d(numPoints).positions

        # Step 2: Convert to Field<Vector2d>
        field = FieldofVector2d("position")
        for x, y in initial:
            field.addValue(Vector2d(x, y))

        print("performing CVT")

        # Step 3: Perform CVT relaxation
        for _ in range(n_iterations):
            vor = VoronoiMesh2d(field)
            cells = vor.getCells()
            assert field.size() == len(cells), "Mismatch in CVT loop"

            for i in range(field.size()):
                centroid = cells[i].centroid
                # Clamp the centroid to stay within [-1, 1] in both x and y
                clamped = Vector2d(
                    max(-1.0, min(1.0, centroid.x)),
                    max(-1.0, min(1.0, centroid.y))
                )
                field.setValue(i, clamped)

        # Step 4: Store result
        self.positions = [[field[i].x, field[i].y] for i in range(field.size())]

class CVTDiskGenerator2d:
    def __init__(self, numPoints, n_iterations=10):
        self.positions = []

        # Step 1: Generate initial seed using Poisson Disk sampling
        from PoissonNodeGenerator import PoissonDisk2d
        initial = PoissonDisk2d(numPoints).positions

        # Step 2: Convert to Field<Vector2d>
        field = FieldofVector2d("position")
        for x, y in initial:
            field.addValue(Vector2d(x, y))

        print("performing CVT in disk")

        # Step 3: Perform CVT relaxation
        for _ in range(n_iterations):
            vor = VoronoiMesh2d(field)
            cells = vor.getCells()
            assert field.size() == len(cells), "Mismatch in CVT loop"

            for i in range(field.size()):
                centroid = cells[i].centroid
                r2 = centroid.mag2
                if r2 > 1.0:
                    centroid *= (1.0 / math.sqrt(r2))  # Project onto unit circle
                field.setValue(i, centroid)

        # Step 4: Store result
        self.positions = [[field[i].x, field[i].y] for i in range(field.size())]
