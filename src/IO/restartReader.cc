// Copyright (C) 2026  Cody Raskin

#ifndef RESTARTREADER_CC
#define RESTARTREADER_CC

#include "restartReader.hh"
#include "restartFormat.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <complex>
#include <stdexcept>

template <typename T>
static void readFieldPayload(std::ifstream& in, NodeList& nodeList, const std::string& name) {
    uint64_t count;
    in.read(reinterpret_cast<char*>(&count), sizeof(count));
    std::vector<T> buffer(count);
    if (count > 0) {
        in.read(reinterpret_cast<char*>(buffer.data()), count * sizeof(T));
    }

    Field<T>* field = nodeList.getField<T>(name);
    if (!field) {
        std::cerr << "RestartReader: field '" << name << "' not present in current NodeList, skipping\n";
        return;
    }
    if (field->size() != count) {
        std::cerr << "RestartReader: size mismatch for field '" << name << "' ("
                  << count << " in file vs " << field->size() << " in NodeList), skipping\n";
        return;
    }
    std::memcpy(field->data(), buffer.data(), count * sizeof(T));
}

// Reads the "id" field payload and reports whether it matches the current
// NodeList's ids exactly -- our sanity check that the restart file lines up
// with the freshly-generated NodeList before we overwrite anything else.
static bool readAndValidateIds(std::ifstream& in, NodeList& nodeList) {
    uint64_t count;
    in.read(reinterpret_cast<char*>(&count), sizeof(count));
    std::vector<int> buffer(count);
    if (count > 0) {
        in.read(reinterpret_cast<char*>(buffer.data()), count * sizeof(int));
    }
    Field<int>* idField = nodeList.getField<int>("id");
    if (!idField || idField->size() != count) return false;
    for (uint64_t i = 0; i < count; ++i) {
        if (idField->getValue(i) != buffer[i]) return false;
    }
    return true;
}

RestartReader::RestartReader(NodeList& nodeList, IntegratorBase& integrator)
    : nodeList(nodeList), integrator(integrator) {}

void
RestartReader::read(const std::string& fileName) {
    std::ifstream in(fileName, std::ios::binary);
    if (!in) {
        throw std::runtime_error("RestartReader: could not open " + fileName + " for reading");
    }

    char magic[4];
    in.read(magic, sizeof(magic));
    if (std::memcmp(magic, RestartFormat::magic, sizeof(magic)) != 0) {
        throw std::runtime_error("RestartReader: " + fileName + " is not a valid Yggdrasil restart file");
    }
    uint32_t version;
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != RestartFormat::version) {
        std::ostringstream msg;
        msg << "RestartReader: " << fileName << " has unsupported restart format version " << version;
        throw std::runtime_error(msg.str());
    }
    int32_t fileDim;
    in.read(reinterpret_cast<char*>(&fileDim), sizeof(fileDim));
    int32_t localDim = nodeList.inferDim();
    if (fileDim != 0 && localDim != 0 && fileDim != localDim) {
        std::ostringstream msg;
        msg << "RestartReader: dimension mismatch in " << fileName
            << " (file appears to be " << fileDim << "D, NodeList appears to be " << localDim << "D)";
        throw std::runtime_error(msg.str());
    }
    uint32_t partitionCount, partitionIndex;
    in.read(reinterpret_cast<char*>(&partitionCount), sizeof(partitionCount));
    in.read(reinterpret_cast<char*>(&partitionIndex), sizeof(partitionIndex));
    uint32_t cycle;
    in.read(reinterpret_cast<char*>(&cycle), sizeof(cycle));
    double time, dt;
    in.read(reinterpret_cast<char*>(&time), sizeof(time));
    in.read(reinterpret_cast<char*>(&dt), sizeof(dt));
    uint64_t numNodes;
    in.read(reinterpret_cast<char*>(&numNodes), sizeof(numNodes));
    if (numNodes != nodeList.size()) {
        std::ostringstream msg;
        msg << "RestartReader: node count mismatch in " << fileName << " ("
            << numNodes << " in file vs " << nodeList.size() << " in NodeList)";
        throw std::runtime_error(msg.str());
    }
    uint64_t numFields;
    in.read(reinterpret_cast<char*>(&numFields), sizeof(numFields));

    for (uint64_t i = 0; i < numFields; ++i) {
        uint32_t nameLength;
        in.read(reinterpret_cast<char*>(&nameLength), sizeof(nameLength));
        std::string name(nameLength, '\0');
        in.read(&name[0], nameLength);
        uint32_t typeTag;
        in.read(reinterpret_cast<char*>(&typeTag), sizeof(typeTag));
        FieldDataType type = static_cast<FieldDataType>(typeTag);

        if (name == "id") {
            if (!readAndValidateIds(in, nodeList)) {
                throw std::runtime_error("RestartReader: 'id' field in " + fileName +
                    " does not match the current NodeList; refusing to restore mismatched state");
            }
            continue;
        }

        switch (type) {
            case FieldDataType::Int:     readFieldPayload<int>(in, nodeList, name); break;
            case FieldDataType::Float:   readFieldPayload<float>(in, nodeList, name); break;
            case FieldDataType::Double:  readFieldPayload<double>(in, nodeList, name); break;
            case FieldDataType::Complex: readFieldPayload<std::complex<double>>(in, nodeList, name); break;
            case FieldDataType::Vector1: readFieldPayload<Lin::Vector<1>>(in, nodeList, name); break;
            case FieldDataType::Vector2: readFieldPayload<Lin::Vector<2>>(in, nodeList, name); break;
            case FieldDataType::Vector3: readFieldPayload<Lin::Vector<3>>(in, nodeList, name); break;
            default:
                std::cerr << "RestartReader: unknown field type tag for '" << name << "', skipping\n";
                break;
        }
    }

    integrator.restoreState(cycle, time, dt);
}

#endif // RESTARTREADER_CC
