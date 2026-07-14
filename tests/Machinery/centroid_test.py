from yggdrasil import *
from Mesh import VoronoiMesh2d
from RandomNodeGenerator import RandomNodeGenerator2d
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Step 1: Generate some random seed points
    numPoints = 10
    seeds = RandomNodeGenerator2d(numPoints).positions

    # Step 2: Create a FieldofVector2d
    field = FieldofVector2d("position")
    for x, y in seeds:
        field.addValue(Vector2d(x, y))

    # Step 3: Create Voronoi mesh from seed points
    vor = VoronoiMesh2d(field)
    cells = vor.getCells()

    # Step 4: Attempt to access centroids
    for i, cell in enumerate(cells):
        try:
            c = cell.centroid
            print(f"Cell {i} centroid = ({c})")
        except Exception as e:
            print(f"Error in centroid for cell {i}: {e}")
