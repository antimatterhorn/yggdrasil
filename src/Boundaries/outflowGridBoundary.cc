// Copyright (C) 2025  Cody Raskin

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
        for (int face = 0; face < 2 * dim; ++face) {
            if (!this->activeFaces[face]) continue;
            const std::vector<int>& ghost = boundaryIds[2 * face];     // ghost cells
            const std::vector<int>& inner = boundaryIds[2 * face + 1]; // first interior cells

            std::vector<int> next;
            next.reserve(inner.size());

            for (int idx : inner) {
                auto [ix, iy, iz] = this->grid->indexToCoordinates(idx);
                int di = 0, dj = 0, dk = 0;

                switch (face) {
                    case 0: di = -1; break;  // -x face
                    case 1: di = +1; break;  // +x face
                    case 2: dj = -1; break;  // -y face
                    case 3: dj = +1; break;  // +y face
                    case 4: dk = -1; break;  // -z face
                    case 5: dk = +1; break;  // +z face
                }

                int nextIdx = this->grid->index(ix + di, iy + dj, iz + dk);
                next.push_back(nextIdx);
            }

            ExtrapolateBoundaryData(field, ghost, inner, next);
        }
    }

    template<typename T>
    void ExtrapolateBoundaryData(Field<T>* field,
                                const std::vector<int>& boundaryIds,
                                const std::vector<int>& copyIds,
                                const std::vector<int>& nextInteriorIds) {
        for (size_t i = 0; i < boundaryIds.size(); ++i) {
            T val = field->getValue(copyIds[i]) * 2.0 - field->getValue(nextInteriorIds[i]);
            field->setValue(boundaryIds[i], val);
        }
    }
public:
    using Vector=Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    OutflowGridBoundary(Mesh::Grid<dim>* grid) : 
        GridBoundary<dim>(grid) {
        if (dim == 1) {
            std::vector<int> leftIds = grid->leftMost();  
            std::vector<int> leftcol;
            leftcol.push_back(1);

            std::vector<int> rightIds = grid->rightMost();
            std::vector<int> rightcol;
            rightcol.push_back(rightIds[0]-1);

            boundaryIds.push_back(leftIds);
            boundaryIds.push_back(leftcol);
            boundaryIds.push_back(rightIds);
            boundaryIds.push_back(rightcol);
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
            
            for (int j=0; j<bottomIds.size(); ++j) {
                leftcol.push_back(grid->index(1,j));
                rightcol.push_back(grid->index(bottomIds.size()-2,j));
            }

            for (int i=0; i<leftIds.size(); ++i) {
                topcol.push_back(grid->index(i,1));
                bottomcol.push_back(grid->index(i,leftIds.size()-2));
            }

            boundaryIds.push_back(leftIds);
            boundaryIds.push_back(leftcol);
            boundaryIds.push_back(rightIds);
            boundaryIds.push_back(rightcol);

            boundaryIds.push_back(topIds);
            boundaryIds.push_back(topcol);
            boundaryIds.push_back(bottomIds);
            boundaryIds.push_back(bottomcol);
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

            boundaryIds.push_back(leftIds);
            boundaryIds.push_back(leftcol);
            boundaryIds.push_back(rightIds);
            boundaryIds.push_back(rightcol);

            boundaryIds.push_back(topIds);
            boundaryIds.push_back(topcol);
            boundaryIds.push_back(bottomIds);
            boundaryIds.push_back(bottomcol);

            boundaryIds.push_back(frontIds);
            boundaryIds.push_back(frontcol);
            boundaryIds.push_back(backIds);
            boundaryIds.push_back(backcol);
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
