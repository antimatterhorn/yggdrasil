// Copyright (C) 2026  Cody Raskin

#ifndef RESTARTREADER_HH
#define RESTARTREADER_HH

#include <string>
#include "../DataBase/nodeList.hh"
#include "../Integrators/integrator.hh"

template <int dim>
class RestartReader {
public:
    RestartReader(NodeList& nodeList, Integrator<dim>& integrator);

    void read(const std::string& fileName);

private:
    NodeList& nodeList;
    Integrator<dim>& integrator;
};

#include "restartReader.cc"

#endif // RESTARTREADER_HH
