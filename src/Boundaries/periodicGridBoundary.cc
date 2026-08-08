// Copyright (C) 2026  Cody Raskin

#include <vector>
#include "gridBoundary.hh"

#define IMPLEMENT_APPLY_INTERFACE(type, TypeName, T) \
    void Apply##TypeName(type* field) override { ApplyThis<T>(field); }

// Base class for Grid Boundary
template <int dim>
class PeriodicGridBoundary : public GridBoundary<dim> {
private:
    std::vector<std::vector<int>> boundaryLists;

    template <typename T>
    void ApplyThis(Field<T>* field) {
        #pragma omp parallel for
        for (int i = 0; i < 2 * dim; ++i) {
            if (!this->activeFaces[i]) continue;
            std::vector<int> b = boundaryLists[2*i];
            std::vector<int> c = boundaryLists[2*i+1];
            CopyBoundaryData<T>(field, b, c);
        }
    }

    template <typename T>
    void CopyBoundaryData(Field<T>* field, const std::vector<int>& boundaryIds, const std::vector<int>& copyIds) {     
        for (int i=0;i<boundaryIds.size();++i) {
            int bi = boundaryIds[i];
            int ci = copyIds[i];
            field->setValue(bi,field->getValue(ci));
        }      
    }
public:
    using Vector      = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;
    using Complex     = std::complex<double>;
    using ComplexField= Field<Complex>;

    PeriodicGridBoundary(Mesh::Grid<dim>* grid) : 
        GridBoundary<dim>(grid) {
        if (dim == 1) {
            std::vector<int> xLow = grid->lowMost(0);
            std::vector<int> xLowInt;
            xLowInt.push_back(1);

            std::vector<int> xHigh = grid->highMost(0);
            std::vector<int> xHighInt;
            xHighInt.push_back(xHigh[0]-1);

            boundaryLists.push_back(xLow);
            boundaryLists.push_back(xHighInt);

            boundaryLists.push_back(xHigh);
            boundaryLists.push_back(xLowInt);
        }
        else if (dim == 2) {
            std::vector<int> xLow  = grid->lowMost(0);
            std::vector<int> xHigh = grid->highMost(0);
            std::vector<int> yLow  = grid->lowMost(1);
            std::vector<int> yHigh = grid->highMost(1);

            std::vector<int> xLowInt, xHighInt, yLowInt, yHighInt;

            int Nx = grid->size_x();
            int Ny = grid->size_y();

            for (int j = 0; j < Ny; ++j) {
                xLowInt.push_back(grid->index(1, j));
                xHighInt.push_back(grid->index(Nx - 2, j));
            }

            for (int i = 0; i < Nx; ++i) {
                yLowInt.push_back(grid->index(i, 1));
                yHighInt.push_back(grid->index(i, Ny - 2));
            }

            boundaryLists.push_back(xLow);
            boundaryLists.push_back(xHighInt);

            boundaryLists.push_back(xHigh);
            boundaryLists.push_back(xLowInt);

            boundaryLists.push_back(yLow);
            boundaryLists.push_back(yHighInt);

            boundaryLists.push_back(yHigh);
            boundaryLists.push_back(yLowInt);
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

            // Iteration order here must match Grid<3>::findBoundaries exactly,
            // since CopyBoundaryData pairs boundaryIds[i] with copyIds[i] by
            // position: j outer / k inner for the x faces.
            for (int j = 0; j < Ny; ++j) {
                for (int k = 0; k < Nz; ++k) {
                    xLowInt.push_back(grid->index(1, j, k));
                    xHighInt.push_back(grid->index(Nx - 2, j, k));
                }
            }

            // i outer / k inner for the y faces.
            for (int i = 0; i < Nx; ++i) {
                for (int k = 0; k < Nz; ++k) {
                    yLowInt.push_back(grid->index(i, 1, k));
                    yHighInt.push_back(grid->index(i, Ny - 2, k));
                }
            }

            // i outer / j inner for the z faces.
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Ny; ++j) {
                    zLowInt.push_back(grid->index(i, j, 1));
                    zHighInt.push_back(grid->index(i, j, Nz - 2));
                }
            }

            boundaryLists.push_back(xLow);
            boundaryLists.push_back(xHighInt);

            boundaryLists.push_back(xHigh);
            boundaryLists.push_back(xLowInt);

            boundaryLists.push_back(yLow);
            boundaryLists.push_back(yHighInt);

            boundaryLists.push_back(yHigh);
            boundaryLists.push_back(yLowInt);

            boundaryLists.push_back(zLow);
            boundaryLists.push_back(zHighInt);

            boundaryLists.push_back(zHigh);
            boundaryLists.push_back(zLowInt);
        }
    }
    
    virtual ~PeriodicGridBoundary() = default;

    // Periodic BCs must always come in axis-aligned pairs (left↔right,
    // top↔bottom, front↔back). Activating one face implicitly activates
    // its partner so the ghost-cell copy is always bidirectional.
    virtual void setFaces(const std::vector<std::string>& faces) override {
        GridBoundary<dim>::setFaces(faces);
        for (int i = 0; i < 2 * dim; i += 2)
            if (this->activeFaces[i] || this->activeFaces[i + 1])
                this->activeFaces[i] = this->activeFaces[i + 1] = true;
    }

    IMPLEMENT_APPLY_INTERFACE(ScalarField, Scalar, double)
    IMPLEMENT_APPLY_INTERFACE(VectorField, Vector, Vector)
    IMPLEMENT_APPLY_INTERFACE(ComplexField, Complex, Complex)

    std::vector<std::vector<int>> GetBounds(Mesh::Grid<dim>* grid) {
        std::vector<std::vector<int>> retVector;
        for (int d = 0; d < dim; ++d) {
            retVector.push_back(grid->lowMost(d));
            retVector.push_back(grid->highMost(d));
        }
        return retVector;
    }
};
