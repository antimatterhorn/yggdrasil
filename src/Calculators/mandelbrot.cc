#include <complex>
#include "../Math/vectorMath.hh"
#include "../DataBase/nodeList.hh"
#include "../Mesh/grid.hh"
#include "../Utilities/progressBar.hh"

class Mandelbrot {
private:
    NodeList* nodeList;
    Mesh::Grid<2>* grid;
public:
    using ScalarField = Field<double>;
    using ComplexField = Field<std::complex<double>>;
    using Vector = Lin::Vector<2>;
    using VectorField = Field<Vector>;

    Mandelbrot(NodeList* nodeList, Mesh::Grid<2>* grid) 
        : nodeList(nodeList),grid(grid) {
        nodeList->template insertField<std::complex<double>>("complexPosition");
        nodeList->template insertField<double>("mandelbrot");
    }

    ~Mandelbrot() {};

    void updatePositions() {
        grid->assignPositions(nodeList);
        ComplexField* complexPosition = nodeList->getField<std::complex<double>>("complexPosition");
        VectorField* position = nodeList->getField<Vector>("position");

        int numNodes = nodeList->size();
        for (int i = 0;i < numNodes;++i) {
            auto& pos = position->getValue(i);
            complexPosition->setValue(i,std::complex<double>(pos[0],pos[1]));
        }
    }

    void compute() {
        updatePositions();
        ComplexField* complexPosition = nodeList->getField<std::complex<double>>("complexPosition");
        ScalarField* mandel = nodeList->getField<double>("mandelbrot");

        int numNodes = nodeList->size();
        int lastPercent = -1;
        #pragma omp parallel for
        for (int i = 0; i < numNodes; ++i) {
            std::complex<double> c = complexPosition->getValue(i);
            std::complex<double> z = 0.0;
            const int maxIterations = 100;
            int iter = 0;
            
            while (std::abs(z) <= 2.0 && iter < maxIterations) {
                z = z*z + c;
                ++iter;
            }

            // If we hit max iterations, assume the point is in the set
            double value = (iter == maxIterations) ? 0.0 : iter + 1 - log(log(abs(z))) / log(2);
            value = std::min(value, (double)maxIterations);
            mandel->setValue(i, value);
            int percent = int(100.0 * i / numNodes);
            if (percent > lastPercent) {
                ProgressBar(double(percent) / 100.0, " ");
                lastPercent = percent;
            }
        }
    }
};