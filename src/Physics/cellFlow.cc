// Copyright (C) 2026  Cody Raskin

#include "physics.hh"
#include "../Trees/kdTree.hh"
#include <cmath>
#include <stdexcept>
#include <algorithm>

// Particle-life-style pairwise forces (CellFlow / Clusters / Ventrella lineage).
// Each node carries a "cellType" species index; forceTable[i*numTypes+j] is the
// (generally asymmetric) interaction coefficient applied to type-i particles from
// type-j neighbors. Force magnitude is a short-range inverse-quadratic repulsion
// minus a quadratic attraction, both in units of the per-type effective radius.

template <int dim>
class CellFlowPhysics : public Physics<dim> {
protected:
    int numTypes;
    std::vector<double> forceTable;
    std::vector<double> radiusByType;
    double baseRadius, repulsion, attraction, k, forceMultiplier, damping;
    double dtmin = 1e30;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    CellFlowPhysics(NodeList* nodeList,
                     PhysicalConstants& constants,
                     int numTypes,
                     std::vector<double> forceTable,
                     std::vector<double> radiusByType,
                     double baseRadius,
                     double repulsion,
                     double attraction,
                     double k,
                     double forceMultiplier,
                     double damping) :
        Physics<dim>(nodeList, constants),
        numTypes(numTypes),
        forceTable(forceTable),
        radiusByType(radiusByType),
        baseRadius(baseRadius),
        repulsion(repulsion),
        attraction(attraction),
        k(k),
        forceMultiplier(forceMultiplier),
        damping(damping) {

        if ((int)forceTable.size() != numTypes*numTypes)
            throw std::runtime_error("CellFlowPhysics: forceTable must have numTypes*numTypes entries");
        if ((int)radiusByType.size() != numTypes)
            throw std::runtime_error("CellFlowPhysics: radiusByType must have numTypes entries");

        this->template EnrollFields<int>({"cellType"});
        // ACCUMULATE: contributes dvdt but does not own or finalize velocity.
        // position/velocity are read from the INTEGRATE package's sub-stage state via initialState.
        this->template EnrollStateFields<Vector>({"velocity"}, FieldPolicy::ACCUMULATE);
    }

    ~CellFlowPhysics() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        int numNodes = nodeList->size();

        VectorField* position = initialState->template getField<Vector>("position");
        VectorField* velocity = initialState->template getField<Vector>("velocity");
        VectorField* dvdt     = deriv.template getField<Vector>("velocity");
        Field<int>*  cellType = nodeList->template getField<int>("cellType");

        double local_dtmin = 1e30;

        KDTree tree(position); // rebuilt each call: positions move every sub-stage
        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i=0; i<numNodes; ++i) {
            int ti = cellType->getValue(i);
            double effRadius = baseRadius * radiusByType[ti];
            Vector xi = position->getValue(i);

            std::vector<int> neighbors = tree.findNearestNeighbors(xi, effRadius);
            Vector accel = Vector::zero();
            for (int j : neighbors) {
                if (j == i) continue;
                int tj = cellType->getValue(j);
                Vector rij = position->getValue(j) - xi;
                double dist = rij.magnitude();
                if (dist < 1e-10) continue;

                double r = dist / effRadius;
                double rep = repulsion / (1.0 + (r*k)*(r*k));
                double att = attraction * r*r;
                double a = forceTable[ti*numTypes + tj];
                double f = a * (rep - att) * forceMultiplier;

                accel += rij.normal() * f;
            }

            Vector vi = velocity->getValue(i);
            dvdt->setValue(i, accel - damping*vi);

            double vmag = vi.magnitude();
            if (vmag > 0.0)
                local_dtmin = std::min(local_dtmin, effRadius/vmag);
        }

        dtmin = local_dtmin;
        this->lastDt = dt;
    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 0.2;
        return timestepCoefficient * dtmin;
    }

    virtual std::string name() const override { return "cellFlow"; }
    virtual std::string description() const override {
        return "CellFlow-style particle-life forces: per-type pairwise attraction/repulsion"; }
};
