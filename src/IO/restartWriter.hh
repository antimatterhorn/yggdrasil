// Copyright (C) 2026  Cody Raskin

#ifndef RESTARTWRITER_HH
#define RESTARTWRITER_HH

#include <string>
#include "../DataBase/nodeList.hh"
#include "../Integrators/integrator.hh"

template <int dim>
class RestartWriter {
public:
    RestartWriter(const NodeList& nodeList, Integrator<dim>& integrator);

    void write(const std::string& fileName);

private:
    const NodeList& nodeList;
    Integrator<dim>& integrator;
};

#include "restartWriter.cc"

#endif // RESTARTWRITER_HH
