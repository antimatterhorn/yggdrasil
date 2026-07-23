// Copyright (C) 2026  Cody Raskin

#ifndef SILOMESHWRITER_HH
#define SILOMESHWRITER_HH

#include <silo.h>
#include <vector>
#include <string>
#include <iostream>
#include "../DataBase/nodeList.hh"
#include "../Math/vectorMath.hh"
#include "../Mesh/grid.hh"
#include "../Mesh/alemesh.hh"

template <int dim>
class SiloMeshWriter {
public:
    SiloMeshWriter(const std::string& baseName, const NodeList& nodeList, const std::vector<std::string>& fieldNames);
    SiloMeshWriter(const std::string& baseName, const NodeList& nodeList, const std::vector<std::string>& fieldNames, Mesh::Grid<dim>* grid);
    SiloMeshWriter(const std::string& baseName, const NodeList& nodeList, const std::vector<std::string>& fieldNames, Mesh::ALEMesh<dim>* aleMesh);

    void write(const std::string& fileName);

private:
    std::string baseName;
    const NodeList& nodeList;
    std::vector<std::string> fieldNames;
    Mesh::Grid<dim>* grid;
    Mesh::ALEMesh<dim>* aleMesh;

    // Zone k in the written ucdmesh corresponds to NodeList/cell index
    // zoneOrder[k] -- Silo's zonelist requires same-shaped zones to be
    // contiguous, so zones get grouped by vertex count and this records the
    // resulting permutation for writeUcdFields to follow.
    std::vector<size_t> zoneOrder;

    void writePointMesh(DBfile* dbfile);
    void writeFields(DBfile* dbfile);

    // Writes the NodeList as zone-centered data on a real Silo quadmesh built
    // from the Grid's cell topology, so VisIt renders actual cells instead of
    // points that need glyphing into squares. Requires nodeList.size() ==
    // grid->size() (one NodeList entry per grid cell); returns false without
    // writing anything if that doesn't hold, so the caller can fall back to
    // writePointMesh/writeFields.
    bool writeQuadMesh(DBfile* dbfile);
    void writeQuadFields(DBfile* dbfile);

    // Writes the NodeList as zone-centered data on a real Silo unstructured
    // (UCD) mesh built from the ALEMesh's actual node positions and cell
    // connectivity -- arbitrary polygon zones, not just quads. Requires
    // nodeList.size() == aleMesh->getConnectivityMap().size() (one NodeList
    // entry per ALEMesh cell); returns false without writing anything if
    // that doesn't hold, so the caller can fall back to writePointMesh.
    bool writeUcdMesh(DBfile* dbfile);
    void writeUcdFields(DBfile* dbfile);
};

#include "siloMeshWriter.cc"

#endif // SILOWRITER_HH