// Copyright (C) 2025  Cody Raskin

#pragma once

#include <vector>
#include "../Math/vectorMath.hh"
#include "basisFunction.hh"
#include <stdexcept>
#include <Eigen/Dense>

namespace Mesh {

enum class ElementType {
    Line,        // dim = 1
    Triangle,    // dim = 2
    Quad,
    Tetrahedron, // dim = 3
    Hexahedron
};

template <int dim>
class Element {
public:
    using Vector = Lin::Vector<dim>;

    Element(ElementType type, const std::vector<size_t>& nodeIndices)
        : mType(type), mNodeIndices(nodeIndices) {}

    virtual ~Element() = default;

    virtual double computeArea(const std::vector<Vector>& positions) const = 0;
    virtual Vector computeCentroid(const std::vector<Vector>& positions) const = 0;
    virtual Eigen::MatrixXd computeStiffnessMatrix(const std::vector<Vector>& positions) const = 0;
    virtual Eigen::MatrixXd computeStructuralStiffnessMatrix(const std::vector<Vector>& positions, const Eigen::MatrixXd& D) const = 0;
    virtual Eigen::VectorXd computeStrain(const std::vector<Vector>& positions, const Eigen::VectorXd& ue) const = 0;

    ElementType type() const { return mType; }
    const std::vector<size_t>& nodeIndices() const { return mNodeIndices; }
    virtual const BasisFunction<dim>& getBasisFunction() const = 0;
    virtual Eigen::VectorXd computeLumpedMassMatrix(const std::vector<Vector>& positions) const = 0;

protected:
    ElementType mType;
    std::vector<size_t> mNodeIndices;
};

class TriangleElement : public Element<2> {
public:
    using Vector = Lin::Vector<2>;

    TriangleElement(const std::vector<size_t>& nodeIndices)
        : Element<2>(ElementType::Triangle, nodeIndices) {}

    double 
    computeArea(const std::vector<Vector>& positions) const override {
        const auto& A = positions[mNodeIndices[0]];
        const auto& B = positions[mNodeIndices[1]];
        const auto& C = positions[mNodeIndices[2]];
        return Lin::triangleArea(A, B, C);  
    }

    Vector 
    computeCentroid(const std::vector<Vector>& positions) const override {
        const auto& A = positions[mNodeIndices[0]];
        const auto& B = positions[mNodeIndices[1]];
        const auto& C = positions[mNodeIndices[2]];
        return Lin::triangleCentroid(A, B, C);
    }

    Eigen::MatrixXd 
    computeStiffnessMatrix(const std::vector<Vector>& positions) const override {
        const auto& A = positions[mNodeIndices[0]];
        const auto& B = positions[mNodeIndices[1]];
        const auto& C = positions[mNodeIndices[2]];
    
        // Build coordinate matrix
        double x1 = A.x(), y1 = A.y();
        double x2 = B.x(), y2 = B.y();
        double x3 = C.x(), y3 = C.y();
    
        // Compute area
        double area = Lin::triangleArea(A, B, C);
        if (area <= 0.0) throw std::runtime_error("Degenerate triangle in stiffness matrix computation");
    
        // Shape function gradients
        Eigen::Matrix<double, 2, 3> Bmat;  // 2 rows (∂/∂x and ∂/∂y), 3 shape functions
    
        Bmat(0, 0) = y2 - y3;
        Bmat(0, 1) = y3 - y1;
        Bmat(0, 2) = y1 - y2;
    
        Bmat(1, 0) = x3 - x2;
        Bmat(1, 1) = x1 - x3;
        Bmat(1, 2) = x2 - x1;
    
        Bmat /= (2.0 * area);
    
        // Compute K = A * (B^T * B); caller scales by conductivity
        Eigen::MatrixXd Ke = area * Bmat.transpose() * Bmat;
    
        return Ke;  // 3x3 matrix
    }

    Eigen::MatrixXd
    computeStructuralStiffnessMatrix(const std::vector<Vector>& positions,
                                     const Eigen::MatrixXd& D) const override {
        const auto& A = positions[mNodeIndices[0]];
        const auto& B = positions[mNodeIndices[1]];
        const auto& C = positions[mNodeIndices[2]];

        double x1 = A.x(), y1 = A.y();
        double x2 = B.x(), y2 = B.y();
        double x3 = C.x(), y3 = C.y();

        double area = Lin::triangleArea(A, B, C);
        if (area <= 0.0) throw std::runtime_error("Degenerate triangle in structural stiffness");

        double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
        double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;

        // Strain-displacement matrix B (3×6): [εxx, εyy, γxy] vs [u1x,u1y, u2x,u2y, u3x,u3y]
        Eigen::Matrix<double, 3, 6> Bmat = Eigen::Matrix<double, 3, 6>::Zero();
        Bmat(0, 0) = b1; Bmat(0, 2) = b2; Bmat(0, 4) = b3;
        Bmat(1, 1) = c1; Bmat(1, 3) = c2; Bmat(1, 5) = c3;
        Bmat(2, 0) = c1; Bmat(2, 1) = b1; Bmat(2, 2) = c2;
        Bmat(2, 3) = b2; Bmat(2, 4) = c3; Bmat(2, 5) = b3;
        Bmat /= (2.0 * area);

        return area * Bmat.transpose() * D * Bmat;  // 6×6
    }

    Eigen::VectorXd
    computeStrain(const std::vector<Vector>& positions,
                  const Eigen::VectorXd& ue) const override {
        const auto& A = positions[mNodeIndices[0]];
        const auto& B = positions[mNodeIndices[1]];
        const auto& C = positions[mNodeIndices[2]];

        double x1 = A.x(), y1 = A.y();
        double x2 = B.x(), y2 = B.y();
        double x3 = C.x(), y3 = C.y();

        double area = Lin::triangleArea(A, B, C);
        if (area <= 0.0) throw std::runtime_error("Degenerate triangle in computeStrain");

        double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
        double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;

        Eigen::Matrix<double, 3, 6> Bmat = Eigen::Matrix<double, 3, 6>::Zero();
        Bmat(0, 0) = b1; Bmat(0, 2) = b2; Bmat(0, 4) = b3;
        Bmat(1, 1) = c1; Bmat(1, 3) = c2; Bmat(1, 5) = c3;
        Bmat(2, 0) = c1; Bmat(2, 1) = b1; Bmat(2, 2) = c2;
        Bmat(2, 3) = b2; Bmat(2, 4) = c3; Bmat(2, 5) = b3;
        Bmat /= (2.0 * area);

        return Bmat * ue;  // [εxx, εyy, γxy]
    }

    const BasisFunction<2>&
    getBasisFunction() const override {
        static TriangleBasisFunction basis;
        return basis;
    }

    Eigen::VectorXd
    computeLumpedMassMatrix(const std::vector<Vector>& positions) const override {
        // this assumes constant density and heat capacity
        double A = computeArea(positions);
        Eigen::VectorXd lumped(3);
        lumped.setConstant(A / 3.0);
        return lumped;
    }
};

class QuadElement : public Element<2> {
public:
    using Vector = Lin::Vector<2>;

    QuadElement(const std::vector<size_t>& nodeIndices)
        : Element<2>(ElementType::Quad, nodeIndices) {
        if (nodeIndices.size() != 4) {
            throw std::invalid_argument("QuadElement must have exactly 4 nodes");
        }
    }

    double 
    computeArea(const std::vector<Vector>& positions) const override {
        const auto& A = positions[mNodeIndices[0]];
        const auto& B = positions[mNodeIndices[1]];
        const auto& C = positions[mNodeIndices[2]];
        const auto& D = positions[mNodeIndices[3]];
        return Lin::quadArea(A, B, C, D);
    }

    Vector 
    computeCentroid(const std::vector<Vector>& positions) const override {
        const auto& A = positions[mNodeIndices[0]];
        const auto& B = positions[mNodeIndices[1]];
        const auto& C = positions[mNodeIndices[2]];
        const auto& D = positions[mNodeIndices[3]];
        return Lin::quadCentroid(A, B, C, D);
    }

    Eigen::MatrixXd 
    computeStiffnessMatrix(const std::vector<Vector>& positions) const override {
        assert(mNodeIndices.size() == 4 && "QuadElement must have 4 nodes");
    
        const auto& A = positions[mNodeIndices[0]];
        const auto& B = positions[mNodeIndices[1]];
        const auto& C = positions[mNodeIndices[2]];
        const auto& D = positions[mNodeIndices[3]];
    
        std::array<Vector, 4> nodeCoords = {A, B, C, D};
    
        // 2x2 Gauss points and weights for reference square [-1, 1] x [-1, 1]
        constexpr double g = 0.5773502691896257;  // 1 / sqrt(3)
        std::array<std::pair<double, double>, 4> gaussPoints = {{
            {-g, -g}, { g, -g}, { g,  g}, { -g,  g}
        }};
        std::array<double, 4> weights = {1.0, 1.0, 1.0, 1.0};
    
        Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(4, 4);
    
        for (int gp = 0; gp < 4; ++gp) {
            double xi = gaussPoints[gp].first;
            double eta = gaussPoints[gp].second;
    
            // Shape function derivatives w.r.t. reference coords (ξ, η)
            Eigen::Matrix<double, 4, 2> dNdXi;
            dNdXi(0, 0) = -0.25 * (1 - eta); dNdXi(0, 1) = -0.25 * (1 - xi);
            dNdXi(1, 0) =  0.25 * (1 - eta); dNdXi(1, 1) = -0.25 * (1 + xi);
            dNdXi(2, 0) =  0.25 * (1 + eta); dNdXi(2, 1) =  0.25 * (1 + xi);
            dNdXi(3, 0) = -0.25 * (1 + eta); dNdXi(3, 1) =  0.25 * (1 - xi);
    
            // Compute Jacobian matrix J
            Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
            for (int i = 0; i < 4; ++i) {
                J(0, 0) += dNdXi(i, 0) * nodeCoords[i].x();
                J(0, 1) += dNdXi(i, 0) * nodeCoords[i].y();
                J(1, 0) += dNdXi(i, 1) * nodeCoords[i].x();
                J(1, 1) += dNdXi(i, 1) * nodeCoords[i].y();
            }
    
            double detJ = J.determinant();
            if (std::abs(detJ) < 1e-15) throw std::runtime_error("Degenerate quad element in stiffness matrix");

            Eigen::Matrix2d invJ = J.inverse();

            // ∇N_i in physical coordinates: dN/dx = dN/dXi * inv(J)
            Eigen::Matrix<double, 4, 2> dNdX;
            for (int i = 0; i < 4; ++i) {
                Eigen::Vector2d gradXi = dNdXi.row(i);
                Eigen::Vector2d gradX = invJ * gradXi;
                dNdX(i, 0) = gradX(0);
                dNdX(i, 1) = gradX(1);
            }

            // Compute local stiffness contribution at this Gauss point
            Eigen::MatrixXd Kgp = Eigen::MatrixXd::Zero(4, 4);
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    double dot = dNdX(i, 0) * dNdX(j, 0) + dNdX(i, 1) * dNdX(j, 1);
                    Kgp(i, j) = dot;
                }
            }

            Ke += Kgp * std::abs(detJ) * weights[gp];
        }
    
        return Ke;
    }

    const BasisFunction<2>& 
    getBasisFunction() const override {
        static QuadBasisFunction basis;
        return basis;
    }

    Eigen::VectorXd 
    computeLumpedMassMatrix(const std::vector<Vector>& positions) const override {
        // this assumes constant density and heat capacity
        std::array<Vector, 4> nodeCoords = {
            positions[mNodeIndices[0]],
            positions[mNodeIndices[1]],
            positions[mNodeIndices[2]],
            positions[mNodeIndices[3]]
        };

        constexpr double g = 0.5773502691896257;
        std::array<std::pair<double, double>, 4> gaussPoints = {{
            {-g, -g}, { g, -g}, { g,  g}, { -g,  g}
        }};
        std::array<double, 4> weights = {1.0, 1.0, 1.0, 1.0};

        Eigen::VectorXd lumped(4);
        lumped.setZero();

        for (int gp = 0; gp < 4; ++gp) {
            double xi = gaussPoints[gp].first;
            double eta = gaussPoints[gp].second;

            // Shape functions at (xi, eta)
            std::array<double, 4> N = {
                0.25 * (1 - xi) * (1 - eta),
                0.25 * (1 + xi) * (1 - eta),
                0.25 * (1 + xi) * (1 + eta),
                0.25 * (1 - xi) * (1 + eta)
            };

            // Compute Jacobian determinant
            Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
            Eigen::Matrix<double, 4, 2> dNdXi;
            dNdXi(0, 0) = -0.25 * (1 - eta); dNdXi(0, 1) = -0.25 * (1 - xi);
            dNdXi(1, 0) =  0.25 * (1 - eta); dNdXi(1, 1) = -0.25 * (1 + xi);
            dNdXi(2, 0) =  0.25 * (1 + eta); dNdXi(2, 1) =  0.25 * (1 + xi);
            dNdXi(3, 0) = -0.25 * (1 + eta); dNdXi(3, 1) =  0.25 * (1 - xi);

            for (int i = 0; i < 4; ++i) {
                J(0, 0) += dNdXi(i, 0) * nodeCoords[i].x();
                J(0, 1) += dNdXi(i, 0) * nodeCoords[i].y();
                J(1, 0) += dNdXi(i, 1) * nodeCoords[i].x();
                J(1, 1) += dNdXi(i, 1) * nodeCoords[i].y();
            }

            double detJ = J.determinant();
            if (std::abs(detJ) < 1e-15) throw std::runtime_error("Degenerate quad element in mass matrix");

            for (int i = 0; i < 4; ++i) {
                lumped(i) += weights[gp] * N[i] * std::abs(detJ);
            }
        }

        return lumped;
    }

    Eigen::MatrixXd
    computeStructuralStiffnessMatrix(const std::vector<Vector>& positions,
                                     const Eigen::MatrixXd& D) const override {
        std::array<Vector, 4> nodeCoords = {
            positions[mNodeIndices[0]], positions[mNodeIndices[1]],
            positions[mNodeIndices[2]], positions[mNodeIndices[3]]
        };

        constexpr double g = 0.5773502691896257;
        std::array<std::pair<double,double>, 4> gps = {{{-g,-g},{g,-g},{g,g},{-g,g}}};

        Eigen::MatrixXd Ke = Eigen::MatrixXd::Zero(8, 8);

        for (int gp = 0; gp < 4; ++gp) {
            double xi = gps[gp].first, eta = gps[gp].second;

            Eigen::Matrix<double, 4, 2> dNdXi;
            dNdXi(0,0) = -0.25*(1-eta); dNdXi(0,1) = -0.25*(1-xi);
            dNdXi(1,0) =  0.25*(1-eta); dNdXi(1,1) = -0.25*(1+xi);
            dNdXi(2,0) =  0.25*(1+eta); dNdXi(2,1) =  0.25*(1+xi);
            dNdXi(3,0) = -0.25*(1+eta); dNdXi(3,1) =  0.25*(1-xi);

            Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
            for (int i = 0; i < 4; ++i) {
                J(0,0) += dNdXi(i,0) * nodeCoords[i].x();
                J(0,1) += dNdXi(i,0) * nodeCoords[i].y();
                J(1,0) += dNdXi(i,1) * nodeCoords[i].x();
                J(1,1) += dNdXi(i,1) * nodeCoords[i].y();
            }

            double detJ = J.determinant();
            if (std::abs(detJ) < 1e-15) throw std::runtime_error("Degenerate quad element in structural stiffness");

            Eigen::Matrix2d invJ = J.inverse();

            // Shape function gradients in physical coords
            Eigen::Matrix<double, 4, 2> dNdX;
            for (int i = 0; i < 4; ++i)
                dNdX.row(i) = (invJ * dNdXi.row(i).transpose()).transpose();

            // Strain-displacement matrix B (3×8): [εxx,εyy,γxy] vs [u1x,u1y,...,u4x,u4y]
            Eigen::Matrix<double, 3, 8> Bmat = Eigen::Matrix<double, 3, 8>::Zero();
            for (int i = 0; i < 4; ++i) {
                Bmat(0, 2*i)   = dNdX(i, 0);
                Bmat(1, 2*i+1) = dNdX(i, 1);
                Bmat(2, 2*i)   = dNdX(i, 1);
                Bmat(2, 2*i+1) = dNdX(i, 0);
            }

            Ke += Bmat.transpose() * D * Bmat * std::abs(detJ);  // weight = 1 for all 4 points
        }

        return Ke;  // 8×8
    }

    Eigen::VectorXd
    computeStrain(const std::vector<Vector>& positions,
                  const Eigen::VectorXd& ue) const override {
        std::array<Vector, 4> nodeCoords = {
            positions[mNodeIndices[0]], positions[mNodeIndices[1]],
            positions[mNodeIndices[2]], positions[mNodeIndices[3]]
        };

        constexpr double g = 0.5773502691896257;
        std::array<std::pair<double,double>, 4> gps = {{{-g,-g},{g,-g},{g,g},{-g,g}}};

        Eigen::Vector3d avgStrain = Eigen::Vector3d::Zero();

        for (int gp = 0; gp < 4; ++gp) {
            double xi = gps[gp].first, eta = gps[gp].second;

            Eigen::Matrix<double, 4, 2> dNdXi;
            dNdXi(0,0) = -0.25*(1-eta); dNdXi(0,1) = -0.25*(1-xi);
            dNdXi(1,0) =  0.25*(1-eta); dNdXi(1,1) = -0.25*(1+xi);
            dNdXi(2,0) =  0.25*(1+eta); dNdXi(2,1) =  0.25*(1+xi);
            dNdXi(3,0) = -0.25*(1+eta); dNdXi(3,1) =  0.25*(1-xi);

            Eigen::Matrix2d J = Eigen::Matrix2d::Zero();
            for (int i = 0; i < 4; ++i) {
                J(0,0) += dNdXi(i,0) * nodeCoords[i].x();
                J(0,1) += dNdXi(i,0) * nodeCoords[i].y();
                J(1,0) += dNdXi(i,1) * nodeCoords[i].x();
                J(1,1) += dNdXi(i,1) * nodeCoords[i].y();
            }

            if (std::abs(J.determinant()) < 1e-15) continue;

            Eigen::Matrix<double, 4, 2> dNdX;
            Eigen::Matrix2d invJ = J.inverse();
            for (int i = 0; i < 4; ++i)
                dNdX.row(i) = (invJ * dNdXi.row(i).transpose()).transpose();

            Eigen::Matrix<double, 3, 8> Bmat = Eigen::Matrix<double, 3, 8>::Zero();
            for (int i = 0; i < 4; ++i) {
                Bmat(0, 2*i)   = dNdX(i,0);
                Bmat(1, 2*i+1) = dNdX(i,1);
                Bmat(2, 2*i)   = dNdX(i,1);
                Bmat(2, 2*i+1) = dNdX(i,0);
            }

            avgStrain += Bmat * ue;
        }

        return avgStrain / 4.0;  // centroid-averaged [εxx, εyy, γxy]
    }
};

template <int dim>
std::shared_ptr<Element<dim>> 
createElement(ElementType type, const std::vector<size_t>& nodeIndices) {
    if constexpr (dim == 2) {
        switch (type) {
            case ElementType::Triangle:
                return std::make_shared<TriangleElement>(nodeIndices);
            case ElementType::Quad:
                return std::make_shared<QuadElement>(nodeIndices);
            default:
                throw std::invalid_argument("Unsupported ElementType for 2D mesh");
        }
    } else if constexpr (dim == 3) {
        // add 3D support here later
        throw std::invalid_argument("3D elements not implemented yet");
    } else {
        throw std::invalid_argument("Unsupported spatial dimension");
    }
}    
}
