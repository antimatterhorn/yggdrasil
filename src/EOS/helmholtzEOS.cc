// Copyright (C) 2026  Cody Raskin

#include <vector>
#include <set>
#include <unordered_map>
#include <tuple>
#include <string>
#include <cmath>
#include <stdexcept>
#include "equationOfState.hh"
#include "eosTable.hh"

class HelmholtzEOS : public EquationOfState {
private:
    std::vector<double> logRhoGrid;
    std::vector<double> logUGrid;

    EOSTable PTable;
    EOSTable CsTable;

    PhysicalConstants& constants;
    std::string tableFile;

    void 
    loadTable(const std::string& filename) {
        FILE* file = std::fopen(filename.c_str(), "r");
        if (!file) {
            throw std::runtime_error("Failed to open EOS table: " + filename);
        }

        char line[512];
        std::vector<std::tuple<double, double, double, double>> rows;
        std::set<double> rawRhos, rawUs;

        // Read file into raw vectors
        while (std::fgets(line, sizeof(line), file)) {
            double logRho, logU, P, cs;
            int count = std::sscanf(line, "%lf %lf %lf %lf", &logRho, &logU, &P, &cs);
            if (count != 4) continue;

            rows.emplace_back(logRho, logU, P, cs);
            rawRhos.insert(logRho);
            rawUs.insert(logU);
        }
        std::fclose(file);

        // Deduplicate and sort axis grids
        logRhoGrid.assign(rawRhos.begin(), rawRhos.end());
        logUGrid.assign(rawUs.begin(), rawUs.end());

        const size_t nRho = logRhoGrid.size();
        const size_t nU   = logUGrid.size();

        std::vector<std::vector<double>> PTableData(nRho, std::vector<double>(nU));
        std::vector<std::vector<double>> CsTableData(nRho, std::vector<double>(nU));

        // Build lookup maps from value to index for fast lookup
        std::unordered_map<double, size_t> rhoIndex, uIndex;
        for (size_t i = 0; i < nRho; ++i) rhoIndex[logRhoGrid[i]] = i;
        for (size_t j = 0; j < nU; ++j)   uIndex[logUGrid[j]] = j;

        // Populate tables
        for (const auto& [logRho, logU, P, cs] : rows) {
            auto itRho = rhoIndex.find(logRho);
            auto itU   = uIndex.find(logU);
            if (itRho == rhoIndex.end() || itU == uIndex.end()) continue;

            size_t i = itRho->second;
            size_t j = itU->second;

            PTableData[i][j]  = P;
            CsTableData[i][j] = cs;
        }

        // Assign to EOSTable members
        PTable  = EOSTable(std::move(PTableData), logRhoGrid, logUGrid);
        CsTable = EOSTable(std::move(CsTableData), logRhoGrid, logUGrid);
    }

    void
    computeHelmholtzApprox(double rho, double u, double& P, double& cs) const {
        const double kB = constants.kB();            // erg/K
        const double mH = constants.protonMass();    // g
        const double mu = 0.6;                       // mean molecular weight
        const double h  = constants.planckConstant();
        const double pi = M_PI;
        const double hb = h/(pi*2.0);
        const double me = constants.electronMass();

        // Convert internal energy per mass to temperature (ideal gas approx)
        double T = (2.0 / 3.0) * (mu * mH / kB) * u;

        // Ideal ion pressure
        double P_ion = (rho / (mu * mH)) * kB * T;

        // Degenerate electron pressure estimate (non-relativistic)
        const double coeff = std::pow(3.0 * pi * pi, 2.0 / 3.0);
        constexpr double mu_e = 2.0;
        double K = (coeff*std::pow(hb,2)) / (5.0 * me);
        double n_e = rho / (mu_e * mH);  // # electrons / cm^3
        double EF = (std::pow(3.0 * pi * pi * n_e, 2.0 / 3.0) * hb * hb) / (2.0 * me);

        // Suppress P_deg if not degenerate
        double P_deg = K * std::pow(n_e, 5.0 / 3.0);

        // Smooth suppression factor
        double ratio = kB*T / EF;
        double suppression = 1.0 / (1.0 + std::pow(ratio, 4.0));
        // Total pressure and sound speed
        P = P_ion + suppression* P_deg;
        double gamma_eff = 5.0 / 3.0;
        cs = std::sqrt(gamma_eff * P / rho);
    }

    // Numerically inverts computeHelmholtzApprox for u at fixed rho: no table is built
    // over (rho, P) since P(rho, u) has no closed-form inverse once the degeneracy
    // suppression term is included.
    double
    invertPressureForU(double rho, double Ptarget) const {
        double uLo = std::pow(10.0, logUGrid.front());
        double uHi = std::pow(10.0, logUGrid.back());

        double PLo, PHi, csDummy;
        computeHelmholtzApprox(rho, uLo, PLo, csDummy);
        computeHelmholtzApprox(rho, uHi, PHi, csDummy);

        if (Ptarget <= PLo) return uLo;
        if (Ptarget >= PHi) return uHi;

        // Bisect in log(u) space for the u that reproduces Ptarget.
        for (int iter = 0; iter < 60; ++iter) {
            double uMid = std::sqrt(uLo * uHi);
            double PMid;
            computeHelmholtzApprox(rho, uMid, PMid, csDummy);
            if (PMid < Ptarget) uLo = uMid;
            else                uHi = uMid;
        }
        return std::sqrt(uLo * uHi);
    }

public:
    HelmholtzEOS(const std::string& tableFile, PhysicalConstants& constants) : 
        EquationOfState(constants), constants(constants), tableFile(tableFile) {
        FILE* fp = std::fopen(tableFile.c_str(), "r");
        if (fp) {
            std::fclose(fp);
            std::cout << "Using table file: " << tableFile << std::endl;
            loadTable(tableFile);
        } else {
            std::cout << "Generating table file: " << tableFile << std::endl;
            generateTable();
        }
    }

    virtual ~HelmholtzEOS();

    virtual void 
    setPressure(Field<double>* pressure, Field<double>* density, Field<double>* internalEnergy) const override {
        for (int i = 0; i < pressure->size(); ++i) {
            double logRho = std::log10(density->getValue(i));
            double logU   = std::log10(internalEnergy->getValue(i));
            double P      = PTable.interpolate(logRho, logU);
            pressure->setValue(i, P);
        }
    }
    virtual void
    setInternalEnergy(Field<double>* internalEnergy, Field<double>* density, Field<double>* pressure) const override {
        for (int i = 0; i < internalEnergy->size(); ++i) {
            double rho = density->getValue(i);
            double P   = pressure->getValue(i);
            double U   = invertPressureForU(rho, P);
            internalEnergy->setValue(i, U);
        }
    }
    virtual void 
    setSoundSpeed(Field<double>* soundSpeed, Field<double>* density, Field<double>* internalEnergy) const override {
        for (int i = 0; i < soundSpeed->size(); ++i) {
            double logRho = std::log10(density->getValue(i));
            double logU   = std::log10(internalEnergy->getValue(i));
            double c      = CsTable.interpolate(logRho, logU);
            soundSpeed->setValue(i, c);
        }
    }
    virtual void 
    setTemperature(Field<double>* temperature, Field<double>* density, Field<double>* internalEnergy) const override {
        const double kB = constants.kB();             // erg/K in code units
        const double mH = constants.protonMass();     // g in code units
        constexpr double mu = 0.6;                    // Mean molecular weight (helium-rich plasma) fix this later

        for (int i = 0; i < temperature->size(); ++i) {
            double u = internalEnergy->getValue(i);
            double T = (2.0 / 3.0) * (mu * mH / kB) * u;
            temperature->setValue(i, T);
        }
    }

    virtual void 
    setInternalEnergyFromTemperature(Field<double>* internalEnergy, Field<double>* density, Field<double>* temperature) const override {
        const double kB = constants.kB();             // erg/K in code units
        const double mH = constants.protonMass();     // g in code units
        constexpr double mu = 0.6;                    // Mean molecular weight

        for (int i = 0; i < internalEnergy->size(); ++i) {
            double T = temperature->getValue(i);
            double u = (3.0 / 2.0) * (kB / (mu * mH)) * T;
            internalEnergy->setValue(i, u);
        }
    }


    virtual void 
    setPressure(double* pressure, double* density, double* internalEnergy) const override {
        double logRho = std::log10(*density);
        double logU   = std::log10(*internalEnergy);
        *pressure     = PTable.interpolate(logRho, logU);
    }

    virtual void
    setInternalEnergy(double* internalEnergy, double* density, double* pressure) const override {
        *internalEnergy = invertPressureForU(*density, *pressure);
    }

    virtual void 
    setSoundSpeed(double* soundSpeed, double* density, double* internalEnergy) const override {
        double logRho = std::log10(*density);
        double logU   = std::log10(*internalEnergy);
        *soundSpeed   = CsTable.interpolate(logRho, logU);
    }

    virtual void 
    setTemperature(double* temperature, double* density, double* internalEnergy) const override {
        const double kB = constants.kB();             // erg/K in code units
        const double mH = constants.protonMass();     // g in code units
        constexpr double mu = 0.6;                    // mean molecular weight

        *temperature = (2.0 / 3.0) * (mu * mH / kB) * (*internalEnergy);
    }

    virtual void 
    setInternalEnergyFromTemperature(double* internalEnergy, double* density, double* temperature) const override {
        const double kB = constants.kB();             // erg/K in code units
        const double mH = constants.protonMass();     // g in code units
        constexpr double mu = 0.6;                    // mean molecular weight

        *internalEnergy = (3.0 / 2.0) * (kB / (mu * mH)) * (*temperature);
    }


    virtual std::string 
    name() const override {
        return "HelmholtzEOS";
    }

    void 
    generateTable() {
        // Define log-spaced rho and u values
        logRhoGrid = Lin::linspace(-2, 10, 300);  // log10(rho) from 1e-2 to 1e10
        logUGrid   = Lin::linspace(-2, 16, 300);  // log10(u) from 1e-2 to 1e16

        std::vector<std::vector<double>> PTableData;
        std::vector<std::vector<double>> CsTableData;

        // Resize tables
        PTableData.resize(logRhoGrid.size(), std::vector<double>(logUGrid.size()));
        CsTableData = PTableData;

        for (size_t i = 0; i < logRhoGrid.size(); ++i) {
            for (size_t j = 0; j < logUGrid.size(); ++j) {
                double rho = std::pow(10.0, logRhoGrid[i]);
                double u   = std::pow(10.0, logUGrid[j]);

                double P, cs;
                computeHelmholtzApprox(rho, u, P, cs);

                PTableData[i][j]  = P;
                CsTableData[i][j] = cs;
            }
        }

        // Write to file using C stdio
        FILE* fp = std::fopen(tableFile.c_str(), "w");
        if (!fp) throw std::runtime_error("Could not open file for writing: " + tableFile);

        for (size_t i = 0; i < logRhoGrid.size(); ++i) {
            for (size_t j = 0; j < logUGrid.size(); ++j) {
                std::fprintf(fp, "%.10e %.10e %.10e %.10e\n",
                            logRhoGrid[i], logUGrid[j],
                            PTableData[i][j], CsTableData[i][j]);
            }
        }

        std::fclose(fp);

        PTable  = EOSTable(std::move(PTableData), logRhoGrid, logUGrid);
        CsTable = EOSTable(std::move(CsTableData), logRhoGrid, logUGrid);
    }

};

HelmholtzEOS::~HelmholtzEOS() = default;