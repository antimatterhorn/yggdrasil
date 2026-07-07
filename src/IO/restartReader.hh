// Copyright (C) 2026  Cody Raskin

#ifndef RESTARTREADER_HH
#define RESTARTREADER_HH

#include <string>
#include "../DataBase/nodeList.hh"
#include "../Integrators/integratorBase.hh"

class RestartReader {
public:
    RestartReader(NodeList& nodeList, IntegratorBase& integrator);

    void read(const std::string& fileName);

private:
    NodeList& nodeList;
    IntegratorBase& integrator;
};

#include "restartReader.cc"

#endif // RESTARTREADER_HH
