// Copyright (C) 2026  Cody Raskin

#include <vector>
#include <unordered_map>
#include "gridBoundary.hh"
#include "../Math/vectorMath.hh"
#include <complex>

#define IMPLEMENT_APPLY_INTERFACE(type, TypeName, T) \
    void Apply##TypeName(type* field) override { ApplyThis<T>(field); }

// Base class for Grid Boundary
template <int dim>
class DirichletGridBoundary : public GridBoundary<dim> {
public:
    using Vector      = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;
    using Complex     = std::complex<double>;
    using ComplexField= Field<Complex>;

private:
    std::vector<int> ids;
    Mesh::Grid<dim>* grid;

    // Saved initial-condition values at Dirichlet cell indices.
    // Populated by ZeroTimeInitialize after the user sets ICs.
    // Restored every step so obstacle cells never become a vacuum.
    std::unordered_map<std::string, std::vector<double>> savedScalars;
    std::unordered_map<std::string, std::vector<Vector>> savedVectors;
    bool initialized = false;

    void
    addIds(std::vector<int> vec) {
        for (int i=0; i<(int)vec.size();i++) {
            ids.push_back(vec[i]);
        }
    }

    template <typename T>
    void ApplyThis(Field<T>* field) {
        if (!initialized) return;
        const std::string name = field->getNameString();
        if constexpr (std::is_same_v<T, double>) {
            auto it = savedScalars.find(name);
            if (it == savedScalars.end()) return;
            const auto& vals = it->second;
            #pragma omp parallel for
            for (int i = 0; i < (int)ids.size(); ++i)
                field->setValue(ids[i], vals[i]);
        } else if constexpr (std::is_same_v<T, Vector>) {
            auto it = savedVectors.find(name);
            if (it == savedVectors.end()) return;
            const auto& vals = it->second;
            #pragma omp parallel for
            for (int i = 0; i < (int)ids.size(); ++i)
                field->setValue(ids[i], vals[i]);
        }
        // Complex fields: not saved, left unchanged
    }

public:
    DirichletGridBoundary(Mesh::Grid<dim>* grid) :
        GridBoundary<dim>(grid),
        grid(grid) {}

    virtual ~DirichletGridBoundary() {}

    // Capture the initial-condition values at all registered Dirichlet cells.
    // Called by Physics::InitializeBoundaries() after the user's IC loop runs.
    virtual void ZeroTimeInitialize(NodeList* nodeList) override {
        std::vector<std::string> names = nodeList->fieldNames();
        for (int fi = 0; fi < (int)names.size(); ++fi) {
            const std::string& name = names[fi];
            FieldBase* fb = nodeList->getFieldByIndex(fi);
            if (auto* sf = dynamic_cast<ScalarField*>(fb)) {
                auto& vals = savedScalars[name];
                vals.resize(ids.size());
                for (int i = 0; i < (int)ids.size(); ++i)
                    vals[i] = sf->getValue(ids[i]);
            } else if (auto* vf = dynamic_cast<VectorField*>(fb)) {
                auto& vals = savedVectors[name];
                vals.resize(ids.size());
                for (int i = 0; i < (int)ids.size(); ++i)
                    vals[i] = vf->getValue(ids[i]);
            }
        }
        initialized = true;
    }

    virtual void
    addBox(Vector p1, Vector p2){
        std::vector<int> local;
        #pragma omp parallel
        {
            std::vector<int> chunk;
            #pragma omp for nowait
            for (int idx = 0; idx < grid->size(); ++idx) {
                Vector thisPos = grid->getPosition(idx);
                bool inside = true;
                for (int i = 0; i < dim; ++i)
                    if (thisPos[i] < p1[i] || thisPos[i] > p2[i]) {
                        inside = false; break;
                    }
                if (inside) chunk.push_back(idx);
            }
            #pragma omp critical
            local.insert(local.end(), chunk.begin(), chunk.end());
        }
        ids.insert(ids.end(), local.begin(), local.end());
    }

    virtual void
    removeBox(Vector p1, Vector p2) {
        std::vector<int> toRemove;
        #pragma omp parallel
        {
            std::vector<int> chunk;
            #pragma omp for nowait
            for (int idx = 0; idx < grid->size(); ++idx) {
                Vector thisPos = grid->getPosition(idx);
                bool inside = true;
                for (int i = 0; i < dim; ++i)
                    if (thisPos[i] < p1[i] || thisPos[i] > p2[i]) {
                        inside = false; break;
                    }
                if (inside) chunk.push_back(idx);
            }
            #pragma omp critical
            toRemove.insert(toRemove.end(), chunk.begin(), chunk.end());
        }
        for (int r : toRemove) {
            auto it = std::find(ids.begin(), ids.end(), r);
            if (it != ids.end()) ids.erase(it);
        }
    }

    virtual void
    addSphere(Vector p, double radius){
        std::vector<int> local;
        #pragma omp parallel
        {
            std::vector<int> chunk;
            #pragma omp for nowait
            for (int idx = 0; idx < grid->size(); ++idx) {
                Vector thisPos = grid->getPosition(idx);
                if ((thisPos - p).mag2() <= radius * radius)
                    chunk.push_back(idx);
            }
            #pragma omp critical
            local.insert(local.end(), chunk.begin(), chunk.end());
        }
        ids.insert(ids.end(), local.begin(), local.end());
    }

    virtual void
    removeSphere(Vector p, double radius){
        std::vector<int> toRemove;
        #pragma omp parallel
        {
            std::vector<int> chunk;
            #pragma omp for nowait
            for (int idx = 0; idx < grid->size(); ++idx) {
                Vector thisPos = grid->getPosition(idx);
                if ((thisPos - p).mag2() <= radius * radius)
                    chunk.push_back(idx);
            }
            #pragma omp critical
            toRemove.insert(toRemove.end(), chunk.begin(), chunk.end());
        }
        for (int r : toRemove) {
            auto it = std::find(ids.begin(), ids.end(), r);
            if (it != ids.end()) ids.erase(it);
        }
    }

    virtual void
    addDomain() {
        for (int d = 0; d < dim; ++d) {
            addIds(grid->lowMost(d));
            addIds(grid->highMost(d));
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
