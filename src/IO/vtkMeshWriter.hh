// Copyright (C) 2026  Cody Raskin

#ifndef VTKMESHWRITER_HH
#define VTKMESHWRITER_HH

#include <vector>
#include <string>
#include <fstream>
#include "../DataBase/nodeList.hh"
#include "../Math/vectorMath.hh"

// Not Grid-aware: writes every NodeList entry as a point, unconditionally.
// For a Grid-backed NodeList (sized NodeList(grid.size())) that includes the
// grid's ghost halo, which shows up here as an extra ring of boundary
// points -- use SiloMeshWriter's grid= quadmesh constructor for a real,
// ghost-free visualization of grid data.
template <int dim>
class VTKMeshWriter {
public:
    VTKMeshWriter(const std::string& baseName, const NodeList* nodeList, const std::vector<std::string>& fieldNames);

    void write(const std::string& fileName);

private:
    std::string baseName;
    const NodeList* nodeList;
    std::vector<std::string> fieldNames;

    void writeVTKHeader(std::ofstream& outFile);
    void writePointCoordinates(std::ofstream& outFile);
    void writePointData(std::ofstream& outFile);
};

#include "vtkMeshWriter.cc"

#endif // VTWMESHWRITER_HH