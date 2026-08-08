# Copyright (C) 2026  Cody Raskin

# Block-structured AMR patch/level bookkeeping, kept in Python (not new C++
# classes) since physics packages have no uniform constructor signature.

from Mesh import Grid2d
from DataBase import NodeList
from LinearAlgebra import Vector2d
from AMR import PatchNeighborBoundary2d


class AMRPatch:
    """One self-contained (Grid, NodeList, Physics, Integrator) tuple."""
    def __init__(self, grid, nodeList, physics, integrator, box, level, nghost):
        self.grid = grid
        self.nodeList = nodeList
        self.physics = physics
        self.integrator = integrator
        self.box = box       # (iLow, jLow, iHigh, jHigh), inclusive interior cells, in this level's own index space
        self.level = level
        self.nghost = nghost


def buildPatch(box, level, dx, dy, nghost, physicsFactory, integratorClass, dtmin, geometry=None):
    """Construct a patch covering `box` (iLow, jLow, iHigh, jHigh, inclusive
    interior cells) at spacing (dx, dy), oversized by `nghost` on every side
    and positioned in world space via Grid.setOrigin."""
    iLow, jLow, iHigh, jHigh = box
    nx = (iHigh - iLow + 1) + 2 * nghost
    ny = (jHigh - jLow + 1) + 2 * nghost

    # World-space position of the oversized grid's lower-left corner.
    worldX = (iLow - nghost) * dx
    worldY = (jLow - nghost) * dy

    grid = Grid2d(nx, ny, dx, dy) if geometry is None else Grid2d(nx, ny, dx, dy, geometry)
    # setOrigin(o) shifts positions by -o, so reaching (worldX, worldY) needs the negated offset.
    grid.setOrigin(Vector2d(-worldX, -worldY))

    nodeList = NodeList(grid.size())
    physics = physicsFactory(nodeList, grid)
    integrator = integratorClass([physics], dtmin=dtmin)

    return AMRPatch(grid, nodeList, physics, integrator, box, level, nghost)


class AMRController:
    def __init__(self, baseGrid, baseNodeList, physicsFactory, integratorClass, dtmin,
                 refinementRatio=2, nghost=1, boundaries=None):
        self.baseGrid = baseGrid
        self.physicsFactory = physicsFactory
        self.integratorClass = integratorClass
        self.dtmin = dtmin
        self.refinementRatio = refinementRatio
        self.nghost = nghost

        basePhysics = physicsFactory(baseNodeList, baseGrid)
        for bc in (boundaries or []):
            basePhysics.addBoundary(bc)
        baseIntegrator = integratorClass([basePhysics], dtmin=dtmin)
        basePatch = AMRPatch(baseGrid, baseNodeList, basePhysics, baseIntegrator,
                              box=(0, 0, baseGrid.nx - 1, baseGrid.ny - 1),
                              level=0, nghost=0)
        self.levels = [[basePatch]]

    def spacing(self, level):
        """(dx, dy) at `level`, level 0 = the base grid's own spacing."""
        ratio = self.refinementRatio ** level
        return self.baseGrid.dx / ratio, self.baseGrid.dy / ratio

    def buildPatch(self, level, box):
        """Construct a new patch at (level, box), register it in self.levels, and return it."""
        dx, dy = self.spacing(level)
        patch = buildPatch(box, level, dx, dy, self.nghost,
                            self.physicsFactory, self.integratorClass, self.dtmin,
                            geometry=self.baseGrid.geometry())
        while len(self.levels) <= level:
            self.levels.append([])
        self.levels[level].append(patch)
        return patch


def _overlapRange(aLo, aHi, bLo, bHi):
    lo, hi = max(aLo, bLo), min(aHi, bHi)
    return (lo, hi) if lo <= hi else None


def _wireFacePairing(nearPatch, farPatch, axis, nghost):
    """If nearPatch is immediately below farPatch along `axis` (0=i, 1=j) with
    overlapping extent on the other axis, register a PatchNeighborBoundary
    filling nearPatch's ghost rim from farPatch's interior. Returns the
    perpendicular overlap range, or None if not adjacent this way."""
    nBox, fBox = nearPatch.box, farPatch.box
    nLo, nHi = nBox[axis], nBox[axis + 2]
    fLo = fBox[axis]
    if nHi + 1 != fLo:
        return None
    other = 1 - axis
    overlap = _overlapRange(nBox[other], nBox[other + 2], fBox[other], fBox[other + 2])
    if overlap is None:
        return None
    perpLo, perpHi = overlap

    nearIds, farIds = [], []
    for p in range(perpLo, perpHi + 1):
        for g in range(nghost):
            nearLocal = [0, 0]
            nearLocal[axis] = nghost + (nHi - nLo) + 1 + g
            nearLocal[other] = nghost + (p - nBox[other])
            farLocal = [0, 0]
            farLocal[axis] = nghost + g
            farLocal[other] = nghost + (p - fBox[other])
            nearIds.append(nearPatch.grid.index(nearLocal[0], nearLocal[1], 0))
            farIds.append(farPatch.grid.index(farLocal[0], farLocal[1], 0))

    nearPatch.physics.addBoundary(PatchNeighborBoundary2d(farPatch.nodeList, nearIds, farIds))
    return (perpLo, perpHi)


def wireSiblingBoundaries(patchA, patchB, nghost):
    """Register PatchNeighborBoundary on both patches for whichever face(s)
    they share (same level, axis-aligned adjacent boxes). Returns False (no
    boundary registered) if they don't share a face."""
    if patchA.level != patchB.level:
        return False
    wired = False
    for axis in (0, 1):
        if _wireFacePairing(patchA, patchB, axis, nghost) is not None:
            wired = True
        if _wireFacePairing(patchB, patchA, axis, nghost) is not None:
            wired = True
    return wired
