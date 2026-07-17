// Copyright (C) 2026  Cody Raskin

#ifndef SILOMESHWRITER_CC
#define SILOMESHWRITER_CC

#include "siloMeshWriter.hh"
#include <map>

template <int dim>
SiloMeshWriter<dim>::SiloMeshWriter(const std::string& baseName, const NodeList& nodeList, const std::vector<std::string>& fieldNames)
    : baseName(baseName), nodeList(nodeList), fieldNames(fieldNames), grid(nullptr), aleMesh(nullptr) {}

template <int dim>
SiloMeshWriter<dim>::SiloMeshWriter(const std::string& baseName, const NodeList& nodeList, const std::vector<std::string>& fieldNames, Mesh::Grid<dim>* grid)
    : baseName(baseName), nodeList(nodeList), fieldNames(fieldNames), grid(grid), aleMesh(nullptr) {}

template <int dim>
SiloMeshWriter<dim>::SiloMeshWriter(const std::string& baseName, const NodeList& nodeList, const std::vector<std::string>& fieldNames, Mesh::ALEMesh<dim>* aleMesh)
    : baseName(baseName), nodeList(nodeList), fieldNames(fieldNames), grid(nullptr), aleMesh(aleMesh) {}

template <int dim>
void
SiloMeshWriter<dim>::write(const std::string& fileName) {
    std::string fullFileName = baseName + fileName;
    DBfile *dbfile = DBCreate(fullFileName.c_str(), DB_CLOBBER, DB_LOCAL, "Simulation data", DB_PDB);

    if (!dbfile) {
        std::cerr << "Error creating Silo file!" << std::endl;
        return;
    }

    if (grid && writeQuadMesh(dbfile)) {
        writeQuadFields(dbfile);
    } else if (aleMesh && writeUcdMesh(dbfile)) {
        writeUcdFields(dbfile);
    } else {
        writePointMesh(dbfile);
        writeFields(dbfile);
    }

    DBClose(dbfile);
}

template <int dim>
void 
SiloMeshWriter<dim>::writePointMesh(DBfile* dbfile) {
    auto positions = nodeList.position<dim>();
    if (!positions) {
        std::cerr << "Position field is not available in the NodeList!" << std::endl;
        return;
    }

    std::vector<float> coords[dim];
    for (int i = 0; i < dim; ++i) {
        coords[i].resize(nodeList.size());
    }

    for (unsigned int idx = 0; idx < nodeList.size(); ++idx) {
        Lin::Vector<dim> pos = positions->getValue(idx);
        for (int i = 0; i < dim; ++i) {
            coords[i][idx] = static_cast<float>(pos[i]);
        }
    }

    float *coords_ptr[dim];
    for (int i = 0; i < dim; ++i) {
        coords_ptr[i] = coords[i].data();
    }

    DBPutPointmesh(dbfile, "pointmesh", dim, coords_ptr, nodeList.size(), DB_FLOAT, nullptr);
}

template <int dim>
void 
SiloMeshWriter<dim>::writeFields(DBfile* dbfile) {
    for (const auto& field_name : fieldNames) {
        if (field_name == "position") continue;  // Skip the position field

        auto field_double = nodeList.getField<double>(field_name);
        if (field_double) {
            std::vector<float> field_data(nodeList.size());
            for (unsigned int i = 0; i < field_double->size(); ++i) {
                field_data[i] = static_cast<float>(field_double->getValue(i));
            }
            DBPutPointvar1(dbfile, field_name.c_str(), "pointmesh", field_data.data(), nodeList.size(), DB_FLOAT, nullptr);
            continue;
        }

        auto field_vector = nodeList.getField<Lin::Vector<dim>>(field_name);
        if (field_vector) {
            std::vector<float> field_data[dim];
            for (int d = 0; d < dim; ++d) {
                field_data[d].resize(nodeList.size());
            }
            for (unsigned int i = 0; i < field_vector->size(); ++i) {
                Lin::Vector<dim> vec = field_vector->getValue(i);
                for (int d = 0; d < dim; ++d) {
                    field_data[d][i] = static_cast<float>(vec[d]);
                }
            }

            float *field_data_ptr[dim];
            for (int d = 0; d < dim; ++d) {
                field_data_ptr[d] = field_data[d].data();
            }

            DBPutPointvar(dbfile, field_name.c_str(), "pointmesh", dim, field_data_ptr, nodeList.size(), DB_FLOAT, nullptr);
        }
    }
}

template <int dim>
bool
SiloMeshWriter<dim>::writeQuadMesh(DBfile* dbfile) {
    if (static_cast<int>(nodeList.size()) != grid->size()) {
        std::cerr << "SiloMeshWriter: NodeList size (" << nodeList.size()
                   << ") does not match grid size (" << grid->size()
                   << "); writing a point mesh instead." << std::endl;
        return false;
    }

    // Grid cells are the NodeList entries (zones); reconstruct the (n+1)
    // corner-node coordinates per axis from the first cell's center, since
    // Grid only stores cell centers and uniform per-axis spacing.
    int ncells[3] = {grid->size_x(), grid->size_y(), grid->size_z()};
    double spacing[3] = {grid->getdx(), grid->getdy(), grid->getdz()};
    Lin::Vector<dim> base = grid->getPosition(0);

    int dims[dim];
    std::vector<float> coords[dim];
    for (int d = 0; d < dim; ++d) {
        dims[d] = ncells[d] + 1;
        coords[d].resize(dims[d]);
        for (int i = 0; i < dims[d]; ++i) {
            coords[d][i] = static_cast<float>(base[d] - 0.5 * spacing[d] + i * spacing[d]);
        }
    }

    float *coords_ptr[dim];
    for (int d = 0; d < dim; ++d) {
        coords_ptr[d] = coords[d].data();
    }

    static const char* axisNames[3] = {"x", "y", "z"};
    const char* coordnames[dim];
    for (int d = 0; d < dim; ++d) {
        coordnames[d] = axisNames[d];
    }

    DBPutQuadmesh(dbfile, "quadmesh", coordnames, coords_ptr, dims, dim, DB_FLOAT, DB_COLLINEAR, nullptr);
    return true;
}

template <int dim>
void
SiloMeshWriter<dim>::writeQuadFields(DBfile* dbfile) {
    int ncells[3] = {grid->size_x(), grid->size_y(), grid->size_z()};
    int zdims[dim];
    for (int d = 0; d < dim; ++d) {
        zdims[d] = ncells[d];
    }

    for (const auto& field_name : fieldNames) {
        if (field_name == "position") continue;  // Skip the position field

        auto field_double = nodeList.getField<double>(field_name);
        if (field_double) {
            std::vector<float> field_data(nodeList.size());
            for (unsigned int i = 0; i < field_double->size(); ++i) {
                field_data[i] = static_cast<float>(field_double->getValue(i));
            }
            DBPutQuadvar1(dbfile, field_name.c_str(), "quadmesh", field_data.data(), zdims, dim, nullptr, 0, DB_FLOAT, DB_ZONECENT, nullptr);
            continue;
        }

        auto field_vector = nodeList.getField<Lin::Vector<dim>>(field_name);
        if (field_vector) {
            std::vector<float> field_data[dim];
            for (int d = 0; d < dim; ++d) {
                field_data[d].resize(nodeList.size());
            }
            for (unsigned int i = 0; i < field_vector->size(); ++i) {
                Lin::Vector<dim> vec = field_vector->getValue(i);
                for (int d = 0; d < dim; ++d) {
                    field_data[d][i] = static_cast<float>(vec[d]);
                }
            }

            float *field_data_ptr[dim];
            static const char* axisSuffix[3] = {"_x", "_y", "_z"};
            std::string componentNames[dim];
            const char* componentNamePtrs[dim];
            for (int d = 0; d < dim; ++d) {
                field_data_ptr[d] = field_data[d].data();
                componentNames[d] = field_name + axisSuffix[d];
                componentNamePtrs[d] = componentNames[d].c_str();
            }

            DBPutQuadvar(dbfile, field_name.c_str(), "quadmesh", dim, componentNamePtrs, field_data_ptr, zdims, dim, nullptr, 0, DB_FLOAT, DB_ZONECENT, nullptr);
        }
    }
}

template <int dim>
bool
SiloMeshWriter<dim>::writeUcdMesh(DBfile* dbfile) {
    const auto& cells = aleMesh->getConnectivityMap();
    const auto& nodes = aleMesh->getNodes();

    if (nodeList.size() != cells.size()) {
        std::cerr << "SiloMeshWriter: NodeList size (" << nodeList.size()
                   << ") does not match ALEMesh cell count (" << cells.size()
                   << "); writing a point mesh instead." << std::endl;
        return false;
    }

    // Node coordinates: the mesh's own node positions, not NodeList's
    // (NodeList is cell-centered here -- one entry per zone).
    std::vector<float> coords[dim];
    for (int d = 0; d < dim; ++d) coords[d].resize(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i)
        for (int d = 0; d < dim; ++d)
            coords[d][i] = static_cast<float>(nodes[i][d]);

    float *coords_ptr[dim];
    for (int d = 0; d < dim; ++d) coords_ptr[d] = coords[d].data();

    static const char* axisNames[3] = {"x", "y", "z"};
    const char* coordnames[dim];
    for (int d = 0; d < dim; ++d) coordnames[d] = axisNames[d];

    // Silo's zonelist requires zones of the same shape (vertex count) to be
    // contiguous, so group cells by vertex count first. std::map keeps the
    // groups in ascending vertex-count order, which is deterministic but
    // otherwise arbitrary -- zoneOrder records which original cell each
    // resulting zone corresponds to, so writeUcdFields can follow the same
    // permutation.
    std::map<size_t, std::vector<size_t>> cellsByShape;
    for (size_t i = 0; i < cells.size(); ++i)
        cellsByShape[cells[i].size()].push_back(i);

    std::vector<int> nodelist;
    std::vector<int> shapetype, shapesize, shapecnt;
    zoneOrder.clear();
    zoneOrder.reserve(cells.size());

    for (const auto& group : cellsByShape) {
        size_t nVerts = group.first;
        const auto& cellIds = group.second;
        if (nVerts < 3) {
            std::cerr << "SiloMeshWriter: skipping " << cellIds.size()
                       << " degenerate cell(s) with < 3 nodes." << std::endl;
            continue;
        }
        for (size_t cellId : cellIds) {
            for (size_t nodeId : cells[cellId])
                nodelist.push_back(static_cast<int>(nodeId));
            zoneOrder.push_back(cellId);
        }
        shapetype.push_back(DB_ZONETYPE_POLYGON);
        shapesize.push_back(static_cast<int>(nVerts));
        shapecnt.push_back(static_cast<int>(cellIds.size()));
    }

    int nzones = static_cast<int>(zoneOrder.size());
    int nshapes = static_cast<int>(shapetype.size());

    DBPutZonelist2(dbfile, "zonelist", nzones, dim,
                    nodelist.data(), static_cast<int>(nodelist.size()), 0,
                    0, 0,
                    shapetype.data(), shapesize.data(), shapecnt.data(), nshapes,
                    nullptr);

    DBPutUcdmesh(dbfile, "ucdmesh", dim, const_cast<char**>(coordnames), coords_ptr,
                 static_cast<int>(nodes.size()), nzones, "zonelist", nullptr, DB_FLOAT, nullptr);

    return true;
}

template <int dim>
void
SiloMeshWriter<dim>::writeUcdFields(DBfile* dbfile) {
    int nzones = static_cast<int>(zoneOrder.size());

    for (const auto& field_name : fieldNames) {
        if (field_name == "position") continue;  // Skip the position field

        auto field_double = nodeList.getField<double>(field_name);
        if (field_double) {
            std::vector<float> field_data(nzones);
            for (int k = 0; k < nzones; ++k)
                field_data[k] = static_cast<float>(field_double->getValue(zoneOrder[k]));
            DBPutUcdvar1(dbfile, field_name.c_str(), "ucdmesh", field_data.data(), nzones, nullptr, 0, DB_FLOAT, DB_ZONECENT, nullptr);
            continue;
        }

        auto field_vector = nodeList.getField<Lin::Vector<dim>>(field_name);
        if (field_vector) {
            std::vector<float> field_data[dim];
            for (int d = 0; d < dim; ++d) field_data[d].resize(nzones);
            for (int k = 0; k < nzones; ++k) {
                Lin::Vector<dim> vec = field_vector->getValue(zoneOrder[k]);
                for (int d = 0; d < dim; ++d)
                    field_data[d][k] = static_cast<float>(vec[d]);
            }

            float *field_data_ptr[dim];
            static const char* axisSuffix[3] = {"_x", "_y", "_z"};
            std::string componentNames[dim];
            const char* componentNamePtrs[dim];
            for (int d = 0; d < dim; ++d) {
                field_data_ptr[d] = field_data[d].data();
                componentNames[d] = field_name + axisSuffix[d];
                componentNamePtrs[d] = componentNames[d].c_str();
            }

            DBPutUcdvar(dbfile, field_name.c_str(), "ucdmesh", dim, componentNamePtrs, field_data_ptr, nzones, nullptr, 0, DB_FLOAT, DB_ZONECENT, nullptr);
        }
    }
}

#endif