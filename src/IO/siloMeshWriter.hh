// Copyright (C) 2025  Cody Raskin

#ifndef SILOMESHWRITER_HH
#define SILOMESHWRITER_HH

#include <silo.h>
#include <vector>
#include <string>
#include <iostream>
#include "../DataBase/nodeList.hh"
#include "../Math/vectorMath.hh"

template <int dim>
class SiloMeshWriter {
public:
    SiloMeshWriter(const std::string& baseName, const NodeList& nodeList, const std::vector<std::string>& fieldNames);

    void write(const std::string& fileName);

private:
    std::string baseName;
    const NodeList& nodeList;
    std::vector<std::string> fieldNames;

    void writePointMesh(DBfile* dbfile);
    void writeFields(DBfile* dbfile);
};

#include "siloMeshWriter.cc"

#endif // SILOWRITER_HH