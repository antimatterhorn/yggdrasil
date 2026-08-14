// Copyright (C) 2026  Cody Raskin

#include <vector>
#include "gridBoundary.hh"

// Base class for Grid Boundary
template <int dim>
class OutflowGridBoundary : public GridBoundary<dim> {
private:
    std::vector<std::vector<int>> boundaryIds;
    std::string derivFieldName;

    template <typename T>
    void ApplyThis(Field<T>* field) {
        // Snapshot so a corner shared by two active faces copies from
        // pre-call values on both, not from the other face's fresh write.
        const std::vector<T> snapshot = field->getValues();

        for (int face = 0; face < 2 * dim; ++face) {
            if (!this->activeFaces[face]) continue;
            const std::vector<int>& ghost = boundaryIds[2 * face];     // ghost cells
            const std::vector<int>& inner = boundaryIds[2 * face + 1]; // first interior cells
            ExtrapolateBoundaryData(field, snapshot, ghost, inner);
        }
    }

    // Zero-gradient outflow: ghost cell copies the first-interior cell.
    template<typename T>
    void ExtrapolateBoundaryData(Field<T>* field,
                                const std::vector<T>& snapshot,
                                const std::vector<int>& boundaryIds,
                                const std::vector<int>& copyIds) {
        assert(copyIds.size() == boundaryIds.size() &&
               "OutflowGridBoundary: ghost and inner cell lists are different sizes");
        for (size_t i = 0; i < boundaryIds.size(); ++i) {
            field->setValue(boundaryIds[i], snapshot[copyIds[i]]);
        }
    }
public:
    using Vector=Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    // boundaryIds[2*face] are the ghost cells; boundaryIds[2*face+1] is the
    // real first/last logical cell they copy from (index 0 / n-1).
    OutflowGridBoundary(Mesh::Grid<dim>* grid) :
        GridBoundary<dim>(grid) {
        if (dim == 1) {
            std::vector<int> xLow = grid->lowMost(0);
            std::vector<int> xLowInt;
            xLowInt.push_back(grid->index(0));

            std::vector<int> xHigh = grid->highMost(0);
            std::vector<int> xHighInt;
            xHighInt.push_back(grid->index(grid->size_x() - 1));

            boundaryIds.push_back(xLow);
            boundaryIds.push_back(xLowInt);
            boundaryIds.push_back(xHigh);
            boundaryIds.push_back(xHighInt);
        }
        else if (dim == 2) {
            std::vector<int> xLow  = grid->lowMost(0);
            std::vector<int> xHigh = grid->highMost(0);
            std::vector<int> yLow  = grid->lowMost(1);
            std::vector<int> yHigh = grid->highMost(1);

            std::vector<int> xLowInt, xHighInt, yLowInt, yHighInt;

            int Nx = grid->size_x();
            int Ny = grid->size_y();

            // x inner cols: one entry per row (ny entries)
            for (int j = 0; j < Ny; ++j) {
                xLowInt.push_back(grid->index(0, j));
                xHighInt.push_back(grid->index(Nx - 1, j));
            }

            // y inner rows: one entry per column (nx entries)
            for (int i = 0; i < Nx; ++i) {
                yLowInt.push_back(grid->index(i, 0));
                yHighInt.push_back(grid->index(i, Ny - 1));
            }

            boundaryIds.push_back(xLow);
            boundaryIds.push_back(xLowInt);
            boundaryIds.push_back(xHigh);
            boundaryIds.push_back(xHighInt);

            boundaryIds.push_back(yLow);
            boundaryIds.push_back(yLowInt);
            boundaryIds.push_back(yHigh);
            boundaryIds.push_back(yHighInt);
        }
        else if (dim == 3) {
            std::vector<int> xLow  = grid->lowMost(0);
            std::vector<int> xHigh = grid->highMost(0);
            std::vector<int> yLow  = grid->lowMost(1);
            std::vector<int> yHigh = grid->highMost(1);
            std::vector<int> zLow  = grid->lowMost(2);
            std::vector<int> zHigh = grid->highMost(2);

            std::vector<int> xLowInt, xHighInt, yLowInt, yHighInt, zLowInt, zHighInt;

            int Nx = grid->size_x();
            int Ny = grid->size_y();
            int Nz = grid->size_z();

            // Iteration order here must match Grid<3>::buildGhostLists exactly,
            // since these pair positionally with the face lists: j outer / k
            // inner for the x faces.
            for (int j = 0; j < Ny; ++j) {
                for (int k = 0; k < Nz; ++k) {
                    xLowInt.push_back(grid->index(0, j, k));
                    xHighInt.push_back(grid->index(Nx - 1, j, k));
                }
            }

            // i outer / k inner for the y faces.
            for (int i = 0; i < Nx; ++i) {
                for (int k = 0; k < Nz; ++k) {
                    yLowInt.push_back(grid->index(i, 0, k));
                    yHighInt.push_back(grid->index(i, Ny - 1, k));
                }
            }

            // i outer / j inner for the z faces.
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    zLowInt.push_back(grid->index(i, j, 0));
                    zHighInt.push_back(grid->index(i, j, Nz - 1));
                }
            }

            boundaryIds.push_back(xLow);
            boundaryIds.push_back(xLowInt);
            boundaryIds.push_back(xHigh);
            boundaryIds.push_back(xHighInt);

            boundaryIds.push_back(yLow);
            boundaryIds.push_back(yLowInt);
            boundaryIds.push_back(yHigh);
            boundaryIds.push_back(yHighInt);

            boundaryIds.push_back(zLow);
            boundaryIds.push_back(zLowInt);
            boundaryIds.push_back(zHigh);
            boundaryIds.push_back(zHighInt);
        }
    }

    OutflowGridBoundary(Mesh::Grid<dim>* grid, std::string derivative) : 
        OutflowGridBoundary<dim>(grid) { 
        derivFieldName = derivative;
        std::cout << "OutflowGridBoundary: " << derivFieldName << std::endl;
    }
    
    virtual ~OutflowGridBoundary() = default;

    virtual void
    ApplyBoundaries(State<dim>* state, NodeList* nodeList) override {
        Mesh::Grid<dim>* grid = this->grid;

        for (int i = 0; i < state->count(); ++i) {
            FieldBase* field = state->getFieldByIndex(i); // Get the field at index i
            if (typeid(*field) == typeid(ScalarField)) {
                ScalarField* doubleField = dynamic_cast<ScalarField*>(field);
                if (doubleField) {
                    ApplyThis<double>(doubleField);
                }
            } else if (typeid(*field) == typeid(VectorField)) {
                VectorField* vectorField = dynamic_cast<VectorField*>(field);
                if (vectorField) {
                    ApplyThis<Vector>(vectorField);
                }
            }
        }

        if (!derivFieldName.empty()) {
            for (int i = 0; i < state->count(); ++i) {
                FieldBase* field = state->getFieldByIndex(i);
                if (field->getNameString() == derivFieldName) {
                    if (typeid(*field) == typeid(ScalarField)) {
                        ScalarField* doubleField = dynamic_cast<ScalarField*>(field);
                        if (doubleField) {
                            ZeroBoundaryDerivative(doubleField, grid);
                        }
                    } else if (typeid(*field) == typeid(VectorField)) {
                        VectorField* vectorField = dynamic_cast<VectorField*>(field);
                        if (vectorField) {
                            ZeroBoundaryDerivative(vectorField, grid);
                        }
                    }
                }
            }
        }
    }

    void ZeroBoundaryDerivative(ScalarField* field, Mesh::Grid<dim>* grid) {
        std::vector<std::vector<int>> bounds = GetBounds(grid);
        for (int d = 0; d<dim;++d) {
            std::vector<int> b1 = bounds[2*d];
            std::vector<int> b2 = bounds[2*d+1];
            for (int i = 0; i<b1.size(); ++i) {
                field->setValue(b1[i],0);
                field->setValue(b2[i],0);
            }
        }
    }

    void ZeroBoundaryDerivative(VectorField* field, Mesh::Grid<dim>* grid) {
        std::vector<std::vector<int>> bounds = GetBounds(grid);
        for (int d = 0; d<dim;++d) {
            std::vector<int> b1 = bounds[2*d];
            std::vector<int> b2 = bounds[2*d+1];
            for (int i = 0; i<b1.size(); ++i) {
                field->setValue(b1[i],Vector());
                field->setValue(b2[i],Vector());
            }
        }
    }

    std::vector<std::vector<int>> GetBounds(Mesh::Grid<dim>* grid) {
        std::vector<std::vector<int>> retVector;
        for (int d = 0; d < dim; ++d) {
            retVector.push_back(grid->lowMost(d));
            retVector.push_back(grid->highMost(d));
        }
        return retVector;
    }
};
