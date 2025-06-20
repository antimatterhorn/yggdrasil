// Copyright (C) 2025  Cody Raskin

#include "physics.hh"
#include "../Mesh/grid.hh"
#include "../IO/importDepthMap.hh"
#include <iostream>

template <int dim>
class WaveEquation : public Physics<dim> {
protected:
    Mesh::Grid<dim>* grid;
    Mesh::Grid<2>* grid2d;
    bool ocean = false;
    double C;
    double dtmin;
    double dxmin = 1e30;
    std::vector<int> insideIds;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    WaveEquation(NodeList* nodeList, PhysicalConstants& constants, Mesh::Grid<dim>* grid, double C) : 
        Physics<dim>(nodeList,constants),
        grid(grid), C(C) {
        VerifyWaveFields();

        grid->assignPositions(nodeList);

        ScalarField* cs = nodeList->getField<double>("soundSpeed");
        for (int i=0; i<nodeList->getNumNodes();++i) cs->setValue(i,C);
        for (int i = 0; i < dim; ++i)
            dxmin = std::min(dxmin, grid->spacing(i));
    }

    WaveEquation(NodeList* nodeList, PhysicalConstants& constants, Mesh::Grid<2>* grid, const std::string& depthMap) : 
        Physics<dim>(nodeList, constants),
        grid2d(grid), C(constants.ESurfaceGrav()), ocean(true) {
        if (dim != 2) {
            std::cerr << "Error: This constructor can only be used with dim = 2" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        VerifyWaveFields();

        grid2d->assignPositions(nodeList);

        this->template EnrollFields<double>({"depth"});
        grid2d->template insertField<double>("depth");
        
        ImportDepthMap map(depthMap);
        map.populateDepthField(grid2d);

        ScalarField* depth      = grid2d->template getField<double>("depth");
        ScalarField* nodeDepth  = nodeList->getField<double>("depth");
        nodeDepth->copyValues(depth);

        ScalarField* cs = nodeList->getField<double>("soundSpeed");
        double maxC = 0;
        #pragma omp parallel for reduction(max:maxC)
        for (int i=0; i<nodeList->getNumNodes();++i) {
            double c = (depth->getValue(i) < 0 ? sqrt(C*std::abs(depth->getValue(i))) : 0);
            cs->setValue(i,c);
            maxC = std::max(c,maxC);
        }
        C = maxC;

        for (int i = 0; i < dim; ++i)
            dxmin = std::min(dxmin, grid2d->spacing(i));
    }

    ~WaveEquation() {}

    virtual void
    ZeroTimeInitialize() override {
        int numNodes = this->nodeList->size();
        for (int i=0; i<numNodes; ++i) {
            if (ocean) {
                if (!grid2d->onBoundary(i))
                    insideIds.push_back(i);
            }
            else {
                if (!grid->onBoundary(i))
                    insideIds.push_back(i);
            }
        }

        this->UpdateState();
        this->InitializeBoundaries();
    }

    void
    VerifyWaveFields() {
        this->template EnrollFields<double>({"phi", "xi", "maxphi", "phisq", "soundSpeed", "waveEnergyDensity"});
        this->template EnrollStateFields<double>({"phi", "xi"});
    }

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {  
        int numNodes = this->nodeList->size();
        
        ScalarField* xi     = initialState->template getField<double>("xi");
        ScalarField* phi    = initialState->template getField<double>("phi");

        ScalarField* DxiDt  = deriv.template getField<double>("xi");
        ScalarField* DphiDt = deriv.template getField<double>("phi");

        ScalarField* cs     = this->nodeList->template getField<double>("soundSpeed");
        ScalarField* e      = this->nodeList->template getField<double>("waveEnergyDensity");

        double local_dtmin = 1e30;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int p = 0; p < insideIds.size(); ++p) {
            int i = insideIds[p];

            double c        = cs->getValue(i);
            double xi_i     = xi->getValue(i);

            double laplace  = grid->laplacian(i, phi);
            auto grad       = grid->gradient(i, phi);
            double grad2    = grad.mag2(); 

            DxiDt->setValue(i, laplace * c * c);
            DphiDt->setValue(i, dt * DxiDt->getValue(i) + xi_i);
            e->setValue(i, 0.5 * (xi_i * xi_i + c * c * grad2));

            local_dtmin = std::min(local_dtmin, 0.2 * dxmin / c);
        }
        dtmin = local_dtmin;
    }

    double
    getCell(int i,int j, std::string fieldName="phi") {
        int idx = (ocean ? grid2d->index(i,j,0) : grid->index(i,j,0));
        ScalarField* phi    = this->nodeList->template getField<double>(fieldName);
        return phi->getValue(idx);
    }

    virtual void
    FinalChecks() override {
        int numNodes = this->nodeList->size();
        
        ScalarField* xi     = this->nodeList->template getField<double>("xi");
        ScalarField* phi    = this->nodeList->template getField<double>("phi");

        ScalarField* mphi   = this->nodeList->template getField<double>("maxphi"); // diagnostic fields for plotting
        ScalarField* phis   = this->nodeList->template getField<double>("phisq");  // diagnostic fields for plotting
        
        for (int i=0; i<numNodes; ++i) {
            mphi->setValue(i,std::max(mphi->getValue(i),phi->getValue(i)*phi->getValue(i)));
            phis->setValue(i,phi->getValue(i)*phi->getValue(i));
        }
    }

    virtual double 
    EstimateTimestep() const override { 
        return dtmin;
    }

    virtual std::string name() const override { return "waveEquation"; }
    virtual std::string description() const override {
        return "Acoustic wave physics package for grids"; }
};


// d^2 phi / dt^2 = c^2 del^2 phi

// IIRC, which solution you get all comes down to initial conditions
// and boundary conditions.

// RHS is second spatial derivative.  Easy.  From Taylor,
// phi_{i+1,j} = phi_ij + dx * phi,x_ij + dx^2/2 phi,,x_ij + ... // where ,,x means d^2/dx^2
// phi_{i-1,j} = phi_ij - dx * phi,x_ij + dx^2/2 phi,,x_ij + ...
// phi_{i,j+1} = phi_ij + dy * phi,y_ij + dy^2/2 phi,,y_ij + ... // where ,,y means d^2/dy^2
// phi_{i,j-1} = phi_ij - dy * phi,y_ij + dy^2/2 phi,,y_ij + ...
// --------- = -----------------------------------
//           = 4 * phi_ij          + 2 * dx^2/2 phi,,x_ij + 2 * dy^2/2 phi,,y_ij + ...
//           = 4 * phi_ij          + dx^2 phi,,x_ij + dy^2 phi,,y_ij + ...
//           = 4 * phi_ij          + dx^2 (phi,,x_ij + phi,,y_ij) + ... // assume dx=dy
//           = 4 * phi_ij          + dx^2 del^2(phi_ij) + ... // assume dx=dy

// del^2 phi_ij = (-4*phi_{i,j} + phi_{i+1,j} + phi_{i-1,j} phi_{i,j+1} + phi_{i,j-1})/dx^2
