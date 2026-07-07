// Copyright (C) 2026  Cody Raskin

#ifndef RESTARTFORMAT_HH
#define RESTARTFORMAT_HH

#include <cstdint>
#include "../DataBase/nodeList.hh"

// Shared binary layout constants for RestartWriter/RestartReader.
//
// File layout:
//   Header:
//     char     magic[4]        "YGRT"
//     uint32_t version
//     int32_t  dim              (best-effort, via NodeList::inferDim(); 0 if unknown)
//     uint32_t partitionCount   (always 1 today; reserved for future MPI decomposition)
//     uint32_t partitionIndex   (always 0 today)
//     uint32_t cycle
//     double   time
//     double   dt
//     uint64_t numNodes
//     uint64_t numFields
//   Then, numFields times:
//     uint32_t nameLength
//     char     name[nameLength]   (not null-terminated)
//     uint32_t typeTag            (FieldDataType)
//     uint64_t count
//     <count elements of the field's element type, raw>
//
// Only numeric field types (int/float/double/complex<double>/Lin::Vector<1,2,3>)
// are written; string and non-numeric fields are skipped since they are not
// evolving simulation state.
//
// Neither RestartWriter nor RestartReader is templated on dim (NodeList isn't
// dim-templated, and IntegratorBase doesn't need to be either), so there's no
// compile-time dim to stamp into the header -- NodeList::inferDim() supplies
// a best-effort value instead (see its definition for the "0 = unknown" case).
namespace RestartFormat {
    constexpr char magic[4] = {'Y', 'G', 'R', 'T'};
    constexpr uint32_t version = 1;
}

#endif // RESTARTFORMAT_HH
