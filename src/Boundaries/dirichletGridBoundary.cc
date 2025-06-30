// Copyright (C) 2025  Cody Raskin

#include <vector>
#include "gridBoundary.hh"
#include "../Math/vectorMath.hh"
#include <complex>

#define IMPLEMENT_APPLY_INTERFACE(type, TypeName, T) \
    void Apply##TypeName(type* field) override { ApplyThis<T>(field); }

// Base class for Grid Boundary
template <int dim>
class DirichletGridBoundary : public GridBoundary<dim> {
private:
    std::vector<int> ids;
    Mesh::Grid<dim>* grid;

    void
    addIds(std::vector<int> vec) {
        for (int i=0; i<vec.size();i++) {
            ids.push_back(vec[i]);
        }
    }

    template <typename T>
    void ApplyThis(Field<T>* field) { 
        for (int i = 0; i < ids.size(); ++i) {
            int k = ids[i];
            if constexpr (std::is_same_v<T, double>) {
                field->setValue(k, 0.0);
            } else {
                field->setValue(k, T());
            }
        }
    }

public:
    using Vector      = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;
    using Complex     = std::complex<double>;
    using ComplexField= Field<Complex>;

    DirichletGridBoundary(Mesh::Grid<dim>* grid) : 
        GridBoundary<dim>(grid),
        grid(grid) {}
    
    virtual ~DirichletGridBoundary() {}

    virtual void
    addBox(Vector p1, Vector p2){
        #pragma omp parallel for
        for (int idx = 0; idx<grid->size(); ++idx) {
            Vector thisPos = grid->getPosition(idx);
            bool inside = true;
            for(int i=0;i<dim;++i)
                if(thisPos[i] < p1[i] || thisPos[i] > p2[i]) {
                    inside = false;
                    break;
                }
            if(inside)
                ids.push_back(idx);
        }
    }

    virtual void 
    removeBox(Vector p1, Vector p2) {
        #pragma omp parallel for
        for (int idx = 0; idx < grid->size(); ++idx) {
            Vector thisPos = grid->getPosition(idx);
            bool inside = true;
            for (int i = 0; i < dim; ++i) {
                if (thisPos[i] < p1[i] || thisPos[i] > p2[i]) {
                    inside = false;
                    break;
                }
            }
            if (inside) {
                // Check if idx is in the ids vector
                auto it = std::find(ids.begin(), ids.end(), idx);
                if (it != ids.end()) {
                    // Element found, remove it
                    #pragma omp critical
                    ids.erase(it);
                }
            }
        }
    }

    virtual void
    addSphere(Vector p, double radius){
        #pragma omp parallel for
        for (int idx = 0; idx<grid->size(); ++idx) {
            Vector thisPos = grid->getPosition(idx);
            bool inside = true;
            if ((thisPos - p).mag2() <= radius*radius)
                ids.push_back(idx);
        }
    }

    virtual void
    removeSphere(Vector p, double radius){
        #pragma omp parallel for
        for (int idx = 0; idx<grid->size(); ++idx) {
            Vector thisPos = grid->getPosition(idx);
            bool inside = true;
            if ((thisPos - p).mag2() <= radius*radius) {
                auto it = std::find(ids.begin(), ids.end(), idx);
                if (it != ids.end()) {
                    // Element found, remove it
                    #pragma omp critical
                    ids.erase(it);
                }
            }              
        }
    }

    virtual void
    addDomain() {
        if (dim == 1) {
            std::vector<int> leftIds = grid->leftMost();  
            std::vector<int> rightIds = grid->rightMost();

            addIds(leftIds);
            addIds(rightIds);
        }
        else if (dim == 2) {
            std::vector<int> leftIds   = grid->leftMost();  
            std::vector<int> rightIds  = grid->rightMost();
            std::vector<int> topIds    = grid->topMost();
            std::vector<int> bottomIds = grid->bottomMost();

            addIds(leftIds);
            addIds(rightIds);
            addIds(topIds);
            addIds(bottomIds);
        }
        else if (dim == 3) {
            std::vector<int> leftIds   = grid->leftMost();  
            std::vector<int> rightIds  = grid->rightMost();
            std::vector<int> topIds    = grid->topMost();
            std::vector<int> bottomIds = grid->bottomMost();
            std::vector<int> frontIds  = grid->frontMost();
            std::vector<int> backIds   = grid->backMost();

            addIds(leftIds);
            addIds(rightIds);
            addIds(topIds);
            addIds(bottomIds);
            addIds(frontIds);
            addIds(backIds);
        }
    }

    virtual std::vector<int> 
    boundaryIds() {
        return ids;
    }

    IMPLEMENT_APPLY_INTERFACE(ScalarField, Scalar, double)
    IMPLEMENT_APPLY_INTERFACE(VectorField, Vector, Vector)
    IMPLEMENT_APPLY_INTERFACE(ComplexField, Complex, Complex)
};
