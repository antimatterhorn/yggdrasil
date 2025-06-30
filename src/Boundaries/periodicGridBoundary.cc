// Copyright (C) 2025  Cody Raskin

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
        for (int i=0; i<2*dim;++i) {
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
            std::vector<int> leftIds = grid->leftMost();  
            std::vector<int> leftcol;
            leftcol.push_back(1);

            std::vector<int> rightIds = grid->rightMost();
            std::vector<int> rightcol;
            rightcol.push_back(rightIds[0]-1);

            boundaryLists.push_back(leftIds);
            boundaryLists.push_back(rightcol); 

            boundaryLists.push_back(rightIds);
            boundaryLists.push_back(leftcol);
            
        }
        else if (dim == 2) {
            std::vector<int> leftIds   = grid->leftMost();  
            std::vector<int> leftcol;

            std::vector<int> rightIds  = grid->rightMost();
            std::vector<int> rightcol;

            std::vector<int> topIds    = grid->topMost();
            std::vector<int> topcol;

            std::vector<int> bottomIds = grid->bottomMost();
            std::vector<int> bottomcol;
            
            for (int j=0; j<rightIds.size(); ++j) {
                leftcol.push_back(grid->index(1,j));
                rightcol.push_back(grid->index(bottomIds.size()-2,j));
            }

            for (int i=0; i<bottomIds.size(); ++i) {
                topcol.push_back(grid->index(i,1));
                bottomcol.push_back(grid->index(i,leftIds.size()-2));
            }

            boundaryLists.push_back(leftIds);
            boundaryLists.push_back(rightcol);

            boundaryLists.push_back(rightIds);
            boundaryLists.push_back(leftcol);

            boundaryLists.push_back(topIds);
            boundaryLists.push_back(bottomcol); 

            boundaryLists.push_back(bottomIds);
            boundaryLists.push_back(topcol);
            
        }
        else if (dim == 3) {
            std::vector<int> leftIds   = grid->leftMost();  
            std::vector<int> leftcol;

            std::vector<int> rightIds  = grid->rightMost();
            std::vector<int> rightcol;

            std::vector<int> topIds    = grid->topMost();
            std::vector<int> topcol;

            std::vector<int> bottomIds = grid->bottomMost();
            std::vector<int> bottomcol;

            std::vector<int> frontIds  = grid->frontMost();
            std::vector<int> frontcol;

            std::vector<int> backIds   = grid->backMost();
            std::vector<int> backcol;

            int Nx = grid->size_x();
            int Ny = grid->size_y();

            for (int k = 0; k < Ny; ++k) {
                for (int j = 0; j < Nx; ++j) {
                    leftcol.push_back(grid->index(1, j, k));
                    rightcol.push_back(grid->index(Nx - 2, j, k));
                }
            }

            for (int k = 0; k < Ny; ++k) {
                for (int i = 0; i < Nx; ++i) {
                    topcol.push_back(grid->index(i, 1, k));
                    bottomcol.push_back(grid->index(i, Ny - 2, k));
                }
            }

            for (int j = 0; j < Nx; ++j) {
                for (int i = 0; i < Ny; ++i) {
                    frontcol.push_back(grid->index(i, j, 1));
                    backcol.push_back(grid->index(i, j, Ny - 2));
                }
            }

            boundaryLists.push_back(leftIds);
            boundaryLists.push_back(rightcol);
            
            boundaryLists.push_back(rightIds);
            boundaryLists.push_back(leftcol);           

            boundaryLists.push_back(topIds);
            boundaryLists.push_back(bottomcol);
          
            boundaryLists.push_back(bottomIds);
            boundaryLists.push_back(topcol);

            boundaryLists.push_back(frontIds);
            boundaryLists.push_back(backcol);
            
            boundaryLists.push_back(backIds);
            boundaryLists.push_back(frontcol);
        }
    }
    
    virtual ~PeriodicGridBoundary() = default;

    IMPLEMENT_APPLY_INTERFACE(ScalarField, Scalar, double)
    IMPLEMENT_APPLY_INTERFACE(VectorField, Vector, Vector)
    IMPLEMENT_APPLY_INTERFACE(ComplexField, Complex, Complex)

    std::vector<std::vector<int>> GetBounds(Mesh::Grid<dim>* grid) {
        std::vector<std::vector<int>> retVector;
        retVector.push_back(grid->leftMost());
        retVector.push_back(grid->rightMost());
        if (dim>1) {
            retVector.push_back(grid->topMost());
            retVector.push_back(grid->bottomMost());
        }
        if (dim>2) {
            retVector.push_back(grid->frontMost());
            retVector.push_back(grid->backMost());
        }
        return retVector;
    }
};
