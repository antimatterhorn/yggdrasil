// Copyright (C) 2026  Cody Raskin

#ifndef RESTARTWRITER_HH
#define RESTARTWRITER_HH

#include <string>
#include "../DataBase/nodeList.hh"
#include "../Integrators/integratorBase.hh"

class RestartWriter {
public:
    RestartWriter(const NodeList& nodeList, IntegratorBase& integrator);

    void write(const std::string& fileName);

private:
    const NodeList& nodeList;
    IntegratorBase& integrator;
};

#include "restartWriter.cc"

#endif // RESTARTWRITER_HH
