from yggdrasil import *
from Mesh import ALEMesh2d

if __name__ == "__main__":
    # 4x4 grid of nodes -> 3x3 grid of unit quad cells, so the center cell
    # is fully surrounded and all 4 of its edges are interior faces.
    mesh = ALEMesh2d()

    def nodeId(i, j):
        return i + j * 4

    def cellId(i, j):
        return i + j * 3

    for j in range(4):
        for i in range(4):
            mesh.addNode(Vector2d(float(i), float(j)))

    for j in range(3):
        for i in range(3):
            mesh.addCell([nodeId(i, j), nodeId(i + 1, j), nodeId(i + 1, j + 1), nodeId(i, j + 1)])

    mesh.computeFaces()
    faces = mesh.getFaces()
    interior = [f for f in faces if not f.isBoundary]
    boundary = [f for f in faces if f.isBoundary]
    print("total faces:", len(faces), " interior:", len(interior), " boundary:", len(boundary))
    # 3x3 cells: 2*3 vertical + 2*3 horizontal interior edges = 12; perimeter = 4*3 = 12 boundary edges.
    assert len(interior) == 12, f"expected 12 interior faces, got {len(interior)}"
    assert len(boundary) == 12, f"expected 12 boundary faces, got {len(boundary)}"

    def facesAround(cell):
        return [f for f in faces if f.leftCell == cell or f.rightCell == cell]

    def divergenceSum(cell):
        sx, sy = 0.0, 0.0
        for f in facesAround(cell):
            # Boundary faces always have the real cell as leftCell (outward
            # normal already); interior faces need the leftCell/rightCell sign flip.
            sign = 1.0 if f.leftCell == cell else -1.0
            sx += sign * f.normal.x * f.area
            sy += sign * f.normal.y * f.area
        return sx, sy

    centerCell = cellId(1, 1)
    touching = facesAround(centerCell)
    print("faces touching center cell:", len(touching))
    assert len(touching) == 4, f"expected 4 faces around the fully-interior center cell, got {len(touching)}"

    # Divergence-theorem sanity check: the outward-oriented normal*area
    # around any closed cell should sum to ~zero -- for the center cell
    # (all 4 faces interior) and for a corner cell (2 interior + 2 boundary
    # faces), verifying the new boundary-face orientation is consistent with
    # the existing interior-face convention.
    sx, sy = divergenceSum(centerCell)
    print("sum of outward normal*area around center cell:", (sx, sy))
    assert abs(sx) < 1e-9 and abs(sy) < 1e-9, "face normals around a closed cell should sum to zero"

    cornerCell = cellId(0, 0)
    cornerFaces = facesAround(cornerCell)
    print("faces touching corner cell (0,0):", len(cornerFaces),
          " boundary:", sum(1 for f in cornerFaces if f.isBoundary))
    assert len(cornerFaces) == 4
    assert sum(1 for f in cornerFaces if f.isBoundary) == 2
    cx, cy = divergenceSum(cornerCell)
    print("sum of outward normal*area around corner cell:", (cx, cy))
    assert abs(cx) < 1e-9 and abs(cy) < 1e-9, "corner cell's interior+boundary faces should also sum to zero"

    centerArea = mesh.cellVolume(centerCell)
    print("center cell area (expect 1.0):", centerArea)
    assert abs(centerArea - 1.0) < 1e-9

    # Move a corner of the center cell and confirm updateGeometry() actually
    # picks up the change (both the cell volume and the moved face).
    movedNode = nodeId(2, 1)
    oldPos = mesh.getNodes()[movedNode]
    mesh.setNodePosition(movedNode, Vector2d(oldPos.x + 0.5, oldPos.y))
    mesh.updateGeometry()

    newArea = mesh.cellVolume(centerCell)
    print("center cell area after moving a corner (expect 1.25):", newArea)
    assert abs(newArea - 1.25) < 1e-9, f"expected updateGeometry() to grow the cell to area 1.25, got {newArea}"

    # The face between the center cell and its right neighbor runs from
    # node(2,1) (the one we moved, now at (2.5,1)) to node(2,2) (unmoved, at
    # (2,2)) -- length sqrt(0.5^2 + 1^2) = sqrt(1.25), not a simple +0.5,
    # since only one of its two endpoints moved and the edge isn't parallel
    # to that motion. getFaces() converts to Python objects on each call, so
    # re-fetch rather than reusing the pre-move `touching` list.
    updatedFaces = mesh.getFaces()
    rightFace = next(f for f in updatedFaces
                      if (f.leftCell, f.rightCell) in [(centerCell, cellId(2, 1)), (cellId(2, 1), centerCell)])
    expectedLength = (0.5 ** 2 + 1.0 ** 2) ** 0.5
    print("right face length after move (expect %.6f):" % expectedLength, rightFace.area)
    assert abs(rightFace.area - expectedLength) < 1e-9

    print("PASS: alemesh_test")
