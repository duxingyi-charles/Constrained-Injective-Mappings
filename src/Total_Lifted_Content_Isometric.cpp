//
// Created by Charles Du on 6/19/21.
//

#include "Total_Lifted_Content_Isometric.h"

//#include <iostream>

void Total_Lifted_Content_Isometric::initialize(const Eigen::MatrixXd &rest_vertices, Eigen::Matrix3Xi faces,
                                                const std::string &form, double alpha) {
    //
    F = faces;
    param_alpha = alpha;
    free_faceI = Eigen::VectorXi::LinSpaced(F.cols(), 0, faces.cols()-1);
//    std::cout << free_faceI << std::endl;

    // compute restD
    compute_squared_edge_Length(rest_vertices,F,restD);

    // compute rest area
    restA.resize(F.cols());
    for (int i = 0; i < F.cols(); ++i) {
        restA(i) = compute_Heron_tri_area(restD(0,i),restD(1,i),restD(2,i));
    }

    double total_rest_area = restA.sum();

    if (form != "harmonic") { // Tutte form
        // in tutte-uniform form,
        // we assume the auxiliary triangles are equilateral triangles
        // whose total measure is the same as the input rest mesh.
        double regular_tri_area = total_rest_area/F.cols();
        double regular_tri_squared_length = 4 * regular_tri_area / sqrt(3);
        restA = Eigen::VectorXd::Constant(F.cols(), regular_tri_area);
        restD = Eigen::MatrixXd::Constant(3,F.cols(),regular_tri_squared_length);
    }

    // compute scaled squared areas of auxiliary triangles: (alpha/2 + alpha^2)*(A_rest)^2
    scaled_squared_restA.resize(F.cols());
    for (int i = 0; i < F.cols(); ++i) {
        scaled_squared_restA(i) = restA(i) * restA(i);
    }
    scaled_squared_restA *= (alpha/2 + alpha*alpha);


    // compute scaled Dirichlet coefficients
    scaled_Dirichlet_coef.resize(3,F.cols());
    for (int i = 0; i < F.cols(); ++i) {
        double d1 = restD(0,i);
        double d2 = restD(1,i);
        double d3 = restD(2,i);
        scaled_Dirichlet_coef(0,i) = d2 + d3 - d1;
        scaled_Dirichlet_coef(1,i) = d3 + d1 - d2;
        scaled_Dirichlet_coef(2,i) = d1 + d2 - d3;
    }
    scaled_Dirichlet_coef *= (alpha/16);

}

double Total_Lifted_Content_Isometric::compute_total_lifted_content_isometric(const Eigen::Matrix2Xd &vertices) const {
    Eigen::VectorXd energyList = Eigen::VectorXd::Zero(F.cols());
    int vDim = 2;
    int simplex_size = 3; //triangle

    for (int i = 0; i < free_faceI.size(); ++i) {
        auto fi = free_faceI(i);
        Eigen::Matrix2Xd vert(vDim, simplex_size);
        vert.col(0) = vertices.col(F(0, fi));
        vert.col(1) = vertices.col(F(1, fi));
        vert.col(2) = vertices.col(F(2, fi));

        energyList(fi) = compute_lifted_TriArea_isometric(vert, scaled_Dirichlet_coef.col(fi),scaled_squared_restA(fi));
    }

    return energyList.sum();
}

double Total_Lifted_Content_Isometric::compute_lifted_TriArea_isometric(const Eigen::Matrix2Xd &vert,
                                                                        const Eigen::Vector3d &Dirichlet_coef,
                                                                        double scaled_squared_rest_area) const {
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto e1 = v2 - v3;
    auto e2 = v3 - v1;
    auto e3 = v1 - v2;

    double area = compute_tri_signed_area(v1,v2,v3);

    return sqrt((1+param_alpha/2)*area*area + Dirichlet_coef(0) * e1.squaredNorm()
    + Dirichlet_coef(1) * e2.squaredNorm() + Dirichlet_coef(2) * e3.squaredNorm()
    + scaled_squared_rest_area
    );
}

double Total_Lifted_Content_Isometric::compute_total_lifted_content_isometric(const Eigen::Matrix2Xd &vertices,
                                                                              Eigen::VectorXd &energyList) const {
    energyList = Eigen::VectorXd::Zero(F.cols());
    int vDim = 2;
    int simplex_size = 3; //triangle

    for (auto i = 0; i < free_faceI.size(); ++i) {
        auto fi = free_faceI(i);
        Eigen::Matrix2Xd vert(vDim, simplex_size);
        vert.col(0) = vertices.col(F(0, fi));
        vert.col(1) = vertices.col(F(1, fi));
        vert.col(2) = vertices.col(F(2, fi));

        energyList(fi) = compute_lifted_TriArea_isometric(
                vert, scaled_Dirichlet_coef.col(fi),scaled_squared_restA(fi));
    }

    return energyList.sum();
}

double
Total_Lifted_Content_Isometric::compute_total_lifted_content_isometric_with_gradient(const Eigen::Matrix2Xd &vertices,
                                                                                     Eigen::Matrix2Xd &grad) const {
    int vDim = 2;
    double energy = 0.0;
    grad = Eigen::Matrix2Xd::Zero(2, vertices.cols());

    for (auto i = 0; i < free_faceI.size(); ++i) {
        auto fi = free_faceI(i);
        int i1, i2, i3;
        i1 = F(0, fi);
        i2 = F(1, fi);
        i3 = F(2, fi);

        Eigen::MatrixXd vert(vDim, 3);
        vert.col(0) = vertices.col(i1);
        vert.col(1) = vertices.col(i2);
        vert.col(2) = vertices.col(i3);

        Eigen::Matrix2Xd g;
        energy += compute_lifted_TriArea_isometric_with_gradient(
                vert,scaled_Dirichlet_coef.col(fi),scaled_squared_restA(fi),g);

        grad.col(i1) += g.col(0);
        grad.col(i2) += g.col(1);
        grad.col(i3) += g.col(2);
    }

    return energy;
}

double Total_Lifted_Content_Isometric::compute_lifted_TriArea_isometric_with_gradient(const Eigen::Matrix2Xd &vert,
                                                                                      const Eigen::Vector3d &Dirichlet_coef,
                                                                                      double scaled_squared_rest_area,
                                                                                      Eigen::Matrix2Xd &grad) const {
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto e1 = v3 - v2;
    auto e2 = v1 - v3;
    auto e3 = v2 - v1;

    double d1 = Dirichlet_coef(0);
    double d2 = Dirichlet_coef(1);
    double d3 = Dirichlet_coef(2);

    double area = compute_tri_signed_area(v1,v2,v3);
    double energy = sqrt((1+param_alpha/2)*area*area + d1 * e1.squaredNorm()
                         + d2 * e2.squaredNorm() + d3 * e3.squaredNorm()
                         + scaled_squared_rest_area);

    auto h1 = rotate_90deg(e1);
    auto h2 = rotate_90deg(e2);
    auto h3 = rotate_90deg(e3);

    double f1 = (0.5+param_alpha/4) * area;

    auto d1e1 = d1 * e1;
    auto d2e2 = d2 * e2;
    auto d3e3 = d3 * e3;

    grad.resize(vert.rows(), vert.cols());
    grad.col(0) = f1 * h1 + (d2e2 - d3e3);
    grad.col(1) = f1 * h2 + (d3e3 - d1e1);
    grad.col(2) = f1 * h3 + (d1e1 - d2e2);

    grad /= energy;

    return energy;
}

