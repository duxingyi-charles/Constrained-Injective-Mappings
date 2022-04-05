//
// Created by Charles Du on 4/5/22.
//

#include "Total_Generalized_Content_3D.h"

#include "geo_util.h"
#include "deformation_gradient_util.h"

void Total_Generalized_Content_3D::initialize(const Eigen::Matrix3Xd &rest_vertices, const Eigen::Matrix4Xi &tets,
                                              const std::string &form, double alpha, double lambda1, double lambda2,
                                              double k) {
    T = tets;
    param_alpha = alpha;
    param_lambda1 = lambda1;
    param_lambda2 = lambda2;
    param_k = k;
    free_tetI = Eigen::VectorXi::LinSpaced(T.cols(), 0, tets.cols()-1);

    // compute rest volumes
    restA.resize(T.cols());
    for (int i = 0; i < T.cols(); ++i) {
        restA(i) = abs(compute_tet_signed_volume(rest_vertices.col(T(0,i)),
                                             rest_vertices.col(T(1,i)),
                                             rest_vertices.col(T(2,i)),
                                             rest_vertices.col(T(3,i))));
    }
    double total_rest_volume = restA.sum();

    if (form != "harmonic") { // Tutte form
        // in tutte-uniform form,
        // we assume the auxiliary tetrahedrons are equilateral tetrahedrons
        // whose total content is the same as the input rest mesh.
        double regular_tet_volume = total_rest_volume/T.cols();
        restA = Eigen::VectorXd::Constant(T.cols(), regular_tet_volume);
        // compute inverse edge matrix and pFpx
        Eigen::Vector3d v1(0,0, 0);
        double regular_tet_edge_length = cbrt(6 * sqrt(2.) * regular_tet_volume);
        Eigen::Vector3d v2(regular_tet_edge_length, 0, 0);
        Eigen::Vector3d v3(regular_tet_edge_length/2, sqrt(3.)*regular_tet_edge_length/2, 0);
        Eigen::Vector3d v4(regular_tet_edge_length/2, sqrt(3.)*regular_tet_edge_length/6,
                           sqrt(6.)*regular_tet_edge_length/3);
        Eigen::Matrix3d edgeMat;
        compute_edge_matrix(v1,v2,v3,v4, edgeMat);
        Eigen::Matrix3d invEdgeMat = edgeMat.inverse();
        Eigen::MatrixXd pFpx;
        compute_flatten_pFpx(invEdgeMat, pFpx);
        rest_invEdgeMat.resize(T.cols(), invEdgeMat);
        pFpx_list.resize(T.cols(), pFpx);
    } else { // harmonic form
        rest_invEdgeMat.reserve(T.cols());
        Eigen::Matrix3d rest_edge_mat;
        for (int i = 0; i < T.cols(); ++i) {
            rest_invEdgeMat.emplace_back();
            compute_edge_matrix(rest_vertices.col(T(0,i)),
                                rest_vertices.col(T(1,i)),
                                rest_vertices.col(T(2,i)),
                                rest_vertices.col(T(3,i)),
                                rest_edge_mat);
            rest_invEdgeMat.back() = rest_edge_mat.inverse();
        }
        // compute derivative of deformation gradient w.r.t. target vertices
        pFpx_list.reserve(T.cols());
        for (int i = 0; i < T.cols(); ++i) {
            pFpx_list.emplace_back();
            compute_flatten_pFpx(rest_invEdgeMat[i], pFpx_list.back());
        }
    }

    // compute squared volumes of auxiliary tets: (A_rest)^2
    squared_restA.resize(T.cols());
    for (int i = 0; i < T.cols(); ++i) {
        squared_restA(i) = restA(i) * restA(i);
    }

    // coefficients
    // 1 + 2 * alpha * lambda2
    energy_constant_1 = 1.0 + 2 * alpha * lambda2;
    // k^(2/3)
    k_2_3 = cbrt(k) * cbrt(k);
    // 2 * alpha * lambda1 * k^(2/3)
    energy_constant_2 = 2.0 * alpha * lambda1 * k_2_3;
    // alpha * lambda1
    energy_constant_3 = alpha * lambda1;
    // 2 * alpha * lambda2 * k^2 + alpha^2
    energy_constant_4 = 2.0 * alpha * lambda2 * k*k + alpha*alpha;
    // 2 * alpha * lambda1
    grad_constant_1 = 2.0 * alpha * lambda1;
}

double Total_Generalized_Content_3D::compute_total_generalized_content(const Eigen::Matrix3Xd &vertices) const {
    Eigen::VectorXd energyList = Eigen::VectorXd::Zero(T.cols());
    int vDim = 3;
    int simplex_size = 4; //tetrahedron

    for (int i = 0; i < free_tetI.size(); ++i) {
        auto ti = free_tetI(i);
        Eigen::Matrix3Xd vert(vDim, simplex_size);
        vert.col(0) = vertices.col(T(0, ti));
        vert.col(1) = vertices.col(T(1, ti));
        vert.col(2) = vertices.col(T(2, ti));
        vert.col(3) = vertices.col(T(3, ti));
        energyList(ti) = compute_generalized_TetVolume(vert, squared_restA(ti),rest_invEdgeMat[ti]);
    }

    return energyList.sum();
}

double Total_Generalized_Content_3D::compute_total_generalized_content(const Eigen::Matrix3Xd &vertices,
                                                                       Eigen::VectorXd &generalized_content_list) const {
    generalized_content_list = Eigen::VectorXd::Zero(T.cols());
    int vDim = 3;
    int simplex_size = 4; //tetrahedron

    for (int i = 0; i < free_tetI.size(); ++i) {
        auto ti = free_tetI(i);
        Eigen::Matrix3Xd vert(vDim, simplex_size);
        vert.col(0) = vertices.col(T(0, ti));
        vert.col(1) = vertices.col(T(1, ti));
        vert.col(2) = vertices.col(T(2, ti));
        vert.col(3) = vertices.col(T(3, ti));
        generalized_content_list(ti) = compute_generalized_TetVolume(vert, squared_restA(ti),rest_invEdgeMat[ti]);
    }

    return generalized_content_list.sum();
}

double
Total_Generalized_Content_3D::compute_generalized_TetVolume(const Eigen::Matrix3Xd &vert, double squared_rest_volume,
                                                            const Eigen::Matrix3d &rest_inverse_EdgeMat) const {
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto v4 = vert.col(3);

    double volume = compute_tet_signed_volume(v1,v2,v3,v4);

    // compute deformation gradient
    Eigen::Matrix3d target_edge_matrix;
    compute_edge_matrix(v1, v2, v3, v4, target_edge_matrix);
    Eigen::Matrix3d deformation_gradient = target_edge_matrix * rest_inverse_EdgeMat;
    // Cauchy-Green Invariants
    Eigen::Matrix3d FTF = deformation_gradient.transpose() * deformation_gradient;
    double I_C = FTF.trace();
    double II_C = FTF.squaredNorm();

    return sqrt(energy_constant_1 * volume * volume +
    squared_rest_volume * (energy_constant_2 * I_C + energy_constant_3 * (I_C*I_C-II_C) + energy_constant_4)
    );
}

double Total_Generalized_Content_3D::compute_total_generalized_content_with_gradient(const Eigen::Matrix3Xd &vertices,
                                                                                     Eigen::Matrix3Xd &grad) const {
    int vDim = 3;
    int simplex_size = 4; // tetrahedron
    double energy = 0.0;
    grad = Eigen::Matrix2Xd::Zero(vDim, vertices.cols());

    for (auto i = 0; i < free_tetI.size(); ++i) {
        auto fi = free_tetI(i);
        int i1, i2, i3, i4;
        i1 = T(0, fi);
        i2 = T(1, fi);
        i3 = T(2, fi);
        i4 = T(3, fi);

        Eigen::MatrixXd vert(vDim, simplex_size);
        vert.col(0) = vertices.col(i1);
        vert.col(1) = vertices.col(i2);
        vert.col(2) = vertices.col(i3);
        vert.col(3) = vertices.col(i4);

        Eigen::Matrix3Xd g;
        energy += compute_generalized_TetVolume_with_gradient(
                vert,squared_restA(fi),rest_invEdgeMat[fi],
                pFpx_list[fi],g);

        grad.col(i1) += g.col(0);
        grad.col(i2) += g.col(1);
        grad.col(i3) += g.col(2);
        grad.col(i4) += g.col(3);
    }

    return energy;
}

double Total_Generalized_Content_3D::compute_generalized_TetVolume_with_gradient(const Eigen::Matrix3Xd &vert,
                                                                                 double squared_rest_volume,
                                                                                 const Eigen::Matrix3d &rest_inverse_EdgeMat,
                                                                                 const Eigen::MatrixXd &pFpx,
                                                                                 Eigen::Matrix3Xd &grad) const {
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto v4 = vert.col(3);

    Eigen::Matrix3Xd dV_dx;
    double volume = compute_tet_signed_volume_with_gradient(v1,v2,v3,v4, dV_dx);

    // compute deformation gradient
    Eigen::Matrix3d target_edge_matrix;
    compute_edge_matrix(v1, v2, v3, v4, target_edge_matrix);
    Eigen::Matrix3d deformation_gradient = target_edge_matrix * rest_inverse_EdgeMat;
    // Cauchy-Green Invariants
    Eigen::Matrix3d FTF = deformation_gradient.transpose() * deformation_gradient;
    double I_C = FTF.trace();
    double II_C = FTF.squaredNorm();

    double energy =  sqrt(energy_constant_1 * volume * volume +
                squared_rest_volume * (energy_constant_2 * I_C + energy_constant_3 * (I_C*I_C-II_C) + energy_constant_4)
    );

    Eigen::Matrix3d pF = (k_2_3 + I_C) * deformation_gradient - deformation_gradient * FTF;
    Eigen::VectorXd vec_pF;
    vectorize(pF, vec_pF);

    grad.resize(vert.rows(), vert.cols());
    grad = energy_constant_1 * volume * dV_dx;
    Eigen::VectorXd px = pFpx.transpose() * (grad_constant_1 * squared_rest_volume * vec_pF);
    grad.col(0) += px.segment<3>(0);
    grad.col(1) += px.segment<3>(4);
    grad.col(2) += px.segment<3>(8);
    grad.col(3) += px.segment<3>(12);
    grad /= energy;

    return energy;
}
