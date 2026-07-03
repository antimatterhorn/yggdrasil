// Copyright (C) 2026  Cody Raskin

#ifndef RESTARTWRITER_CC
#define RESTARTWRITER_CC

#include "restartWriter.hh"
#include "restartFormat.hh"
#include <fstream>
#include <iostream>
#include <complex>

template <typename T>
static void writeFieldPayload(std::ofstream& out, FieldBase* fieldPtr) {
    Field<T>* field = static_cast<Field<T>*>(fieldPtr);
    uint64_t count = field->size();
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));
    out.write(reinterpret_cast<const char*>(field->data()), count * sizeof(T));
}

RestartWriter::RestartWriter(const NodeList& nodeList, IntegratorBase& integrator)
    : nodeList(nodeList), integrator(integrator) {}

void
RestartWriter::write(const std::string& fileName) {
    std::ofstream out(fileName, std::ios::binary);
    if (!out) {
        std::cerr << "RestartWriter: could not open " << fileName << " for writing\n";
        return;
    }

    std::vector<FieldBase*> checkpointable;
    for (size_t i = 0; i < nodeList.getFieldCount(); ++i) {
        FieldBase* f = nodeList.getFieldByIndex(i);
        if (!f || !f->hasName()) continue;
        if (f->type() == FieldDataType::Unknown || f->type() == FieldDataType::String) {
            std::cerr << "RestartWriter: skipping field '" << f->getNameString()
                      << "' (not a checkpointable type)\n";
            continue;
        }
        checkpointable.push_back(f);
    }

    out.write(RestartFormat::magic, sizeof(RestartFormat::magic));
    uint32_t version = RestartFormat::version;
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    int32_t dimVal = nodeList.inferDim();
    out.write(reinterpret_cast<const char*>(&dimVal), sizeof(dimVal));
    uint32_t partitionCount = 1, partitionIndex = 0;
    out.write(reinterpret_cast<const char*>(&partitionCount), sizeof(partitionCount));
    out.write(reinterpret_cast<const char*>(&partitionIndex), sizeof(partitionIndex));
    uint32_t cycle = integrator.Cycle();
    out.write(reinterpret_cast<const char*>(&cycle), sizeof(cycle));
    double time = integrator.Time();
    out.write(reinterpret_cast<const char*>(&time), sizeof(time));
    double dt = integrator.Dt();
    out.write(reinterpret_cast<const char*>(&dt), sizeof(dt));
    uint64_t numNodes = nodeList.size();
    out.write(reinterpret_cast<const char*>(&numNodes), sizeof(numNodes));
    uint64_t numFields = checkpointable.size();
    out.write(reinterpret_cast<const char*>(&numFields), sizeof(numFields));

    for (FieldBase* f : checkpointable) {
        std::string name = f->getNameString();
        uint32_t nameLength = static_cast<uint32_t>(name.size());
        out.write(reinterpret_cast<const char*>(&nameLength), sizeof(nameLength));
        out.write(name.data(), nameLength);
        uint32_t typeTag = static_cast<uint32_t>(f->type());
        out.write(reinterpret_cast<const char*>(&typeTag), sizeof(typeTag));

        switch (f->type()) {
            case FieldDataType::Int:     writeFieldPayload<int>(out, f); break;
            case FieldDataType::Float:   writeFieldPayload<float>(out, f); break;
            case FieldDataType::Double:  writeFieldPayload<double>(out, f); break;
            case FieldDataType::Complex: writeFieldPayload<std::complex<double>>(out, f); break;
            case FieldDataType::Vector1: writeFieldPayload<Lin::Vector<1>>(out, f); break;
            case FieldDataType::Vector2: writeFieldPayload<Lin::Vector<2>>(out, f); break;
            case FieldDataType::Vector3: writeFieldPayload<Lin::Vector<3>>(out, f); break;
            default: break; // unreachable: filtered out above
        }
    }
}

#endif // RESTARTWRITER_CC
