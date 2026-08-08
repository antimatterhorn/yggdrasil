# Copyright (C) 2026  Cody Raskin

# Block-structured AMR patch/level bookkeeping, kept in Python (not new C++
# classes) since physics packages have no uniform constructor signature.

from Mesh import Grid2d
from DataBase import NodeList
from LinearAlgebra import Vector2d
from AMR import (PatchNeighborBoundary2d, CoarseFineBoundary2d, RestrictionOperator2d,
                  FluxRegister2d, Refluxer2d)


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


def wireCoarseFineBoundary(finePatch, coarsePatch, refinementRatio):
    """Register a CoarseFineBoundary on finePatch filling each of its ghost
    cells from the coarse cell of coarsePatch containing it (piecewise-constant
    injection). Ghost cells with no covering coarse cell -- e.g. on a face that
    is really a domain boundary -- are skipped. Returns the number of cells
    paired."""
    iLow, jLow, iHigh, jHigh = finePatch.box
    nghost = finePatch.nghost
    interiorNx = iHigh - iLow + 1
    interiorNy = jHigh - jLow + 1

    cLowI, cLowJ = coarsePatch.box[0], coarsePatch.box[1]
    cGhost = coarsePatch.nghost

    ghostIds, coarseIds = [], []
    for lj in range(interiorNy + 2 * nghost):
        for li in range(interiorNx + 2 * nghost):
            if nghost <= li < nghost + interiorNx and nghost <= lj < nghost + interiorNy:
                continue
            # Floor division is what we want here: a negative fine index still lands in the coarse cell covering it.
            cLocalI = (iLow + li - nghost) // refinementRatio - cLowI + cGhost
            cLocalJ = (jLow + lj - nghost) // refinementRatio - cLowJ + cGhost
            if not (0 <= cLocalI < coarsePatch.grid.nx and 0 <= cLocalJ < coarsePatch.grid.ny):
                continue
            ghostIds.append(finePatch.grid.index(li, lj, 0))
            coarseIds.append(coarsePatch.grid.index(cLocalI, cLocalJ, 0))

    finePatch.physics.addBoundary(CoarseFineBoundary2d(coarsePatch.nodeList, ghostIds, coarseIds))
    return len(ghostIds)


def buildRestriction(finePatch, coarsePatch, refinementRatio, eos):
    """Build a RestrictionOperator averaging finePatch's interior cells onto
    the coarse cells they cover. Only coarse cells covered in full by fine
    interior cells are included -- a partially covered cell (possible when the
    fine box isn't aligned to the coarse grid) would otherwise be overwritten
    using only part of its volume."""
    iLow, jLow, iHigh, jHigh = finePatch.box
    nghost = finePatch.nghost
    cLowI, cLowJ = coarsePatch.box[0], coarsePatch.box[1]
    cGhost = coarsePatch.nghost

    perCoarse = {}
    for lj in range(nghost, nghost + (jHigh - jLow + 1)):
        for li in range(nghost, nghost + (iHigh - iLow + 1)):
            cLocalI = (iLow + li - nghost) // refinementRatio - cLowI + cGhost
            cLocalJ = (jLow + lj - nghost) // refinementRatio - cLowJ + cGhost
            if not (0 <= cLocalI < coarsePatch.grid.nx and 0 <= cLocalJ < coarsePatch.grid.ny):
                continue
            key = coarsePatch.grid.index(cLocalI, cLocalJ, 0)
            perCoarse.setdefault(key, []).append(finePatch.grid.index(li, lj, 0))

    fullyCovered = refinementRatio ** 2
    fineIds, coarseIds = [], []
    for cIdx, fIdxs in perCoarse.items():
        if len(fIdxs) != fullyCovered:
            continue
        for f in fIdxs:
            fineIds.append(f)
            coarseIds.append(cIdx)

    return RestrictionOperator2d(finePatch.nodeList, coarsePatch.nodeList,
                                  finePatch.grid, coarsePatch.grid, eos,
                                  fineIds, coarseIds)


class Reflux:
    """Flux registers for a coarse-fine interface plus the correction that uses
    them. Holds both registers so they outlive the raw pointers handed to the
    Refluxer and to each solver."""
    def __init__(self, refluxer, coarseRegister, fineRegister):
        self.refluxer = refluxer
        self.coarseRegister = coarseRegister
        self.fineRegister = fineRegister

    def apply(self, dt):
        self.refluxer.apply(dt)


def buildReflux(finePatch, coarsePatch, refinementRatio, eos):
    """Register the faces along finePatch's boundary and the coarse faces they
    cover, attach a flux register to each solver, and return the Reflux that
    corrects the coarse cells once both levels have advanced."""
    r = refinementRatio
    box = finePatch.box
    for axis in (0, 1):
        assert box[axis] % r == 0 and (box[axis + 2] + 1) % r == 0, \
            "reflux needs the fine box aligned to coarse cell boundaries"

    nghost = finePatch.nghost
    lowLevel = [box[0], box[1]]
    interior = [box[2] - box[0] + 1, box[3] - box[1] + 1]
    cLow = [coarsePatch.box[0], coarsePatch.box[1]]
    cGhost = coarsePatch.nghost

    coarseRegister = FluxRegister2d(coarsePatch.grid.size())
    fineRegister = FluxRegister2d(finePatch.grid.size())
    coarseSlots, fineSlots = [], []

    for axis in (0, 1):
        other = 1 - axis
        covLow = lowLevel[axis] // r
        covHigh = (lowLevel[axis] + interior[axis] - 1) // r
        perpLow = lowLevel[other] // r
        perpHigh = (lowLevel[other] + interior[other] - 1) // r

        # fineOnPlus: the coarse cell just below the patch, whose +axis face borders it.
        for fineOnPlus in (True, False):
            cAlong = covLow - 1 if fineOnPlus else covHigh + 1
            fineEdge = nghost if fineOnPlus else nghost + interior[axis] - 1
            for cPerp in range(perpLow, perpHigh + 1):
                cLocal = [0, 0]
                cLocal[axis] = cAlong - cLow[axis] + cGhost
                cLocal[other] = cPerp - cLow[other] + cGhost
                if not (0 <= cLocal[0] < coarsePatch.grid.nx and 0 <= cLocal[1] < coarsePatch.grid.ny):
                    continue
                cSlot = coarseRegister.registerFace(coarsePatch.grid.index(cLocal[0], cLocal[1], 0),
                                                     axis, fineOnPlus)
                for k in range(r):
                    fLocal = [0, 0]
                    fLocal[axis] = fineEdge
                    fLocal[other] = cPerp * r + k - lowLevel[other] + nghost
                    fSlot = fineRegister.registerFace(finePatch.grid.index(fLocal[0], fLocal[1], 0),
                                                       axis, not fineOnPlus)
                    coarseSlots.append(cSlot)
                    fineSlots.append(fSlot)

    coarsePatch.physics.attachFluxObserver(coarseRegister)
    finePatch.physics.attachFluxObserver(fineRegister)
    refluxer = Refluxer2d(coarsePatch.nodeList, coarsePatch.grid, eos,
                           coarseRegister, fineRegister, coarseSlots, fineSlots)
    return Reflux(refluxer, coarseRegister, fineRegister)


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
