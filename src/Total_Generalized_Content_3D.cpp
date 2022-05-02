//
// Created by Charles Du on 4/5/22.
//

#include "Total_Generalized_Content_3D.h"

#include <Eigen/Dense> // for inverse()
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

#pragma omp parallel
#pragma omp for
    for (int i = 0; i < free_tetI.size(); ++i) {
        auto ti = free_tetI(i);
        energyList(ti) = compute_generalized_TetVolume(vertices.col(T(0, ti)),
                                                       vertices.col(T(1, ti)),
                                                       vertices.col(T(2, ti)),
                                                       vertices.col(T(3, ti)),
                                                       squared_restA(ti),rest_invEdgeMat[ti]);
    }

    return energyList.sum();
}

double Total_Generalized_Content_3D::compute_total_generalized_negative_content(const Eigen::Matrix3Xd &vertices) const {
    Eigen::VectorXd energyList = Eigen::VectorXd::Zero(T.cols());
    int vDim = 3;
    int simplex_size = 4; //tetrahedron

#pragma omp parallel
#pragma omp for
    for (int i = 0; i < free_tetI.size(); ++i) {
        auto ti = free_tetI(i);
        energyList(ti) = compute_generalized_negative_TetVolume(vertices.col(T(0, ti)),
                                                       vertices.col(T(1, ti)),
                                                       vertices.col(T(2, ti)),
                                                       vertices.col(T(3, ti)),
                                                       squared_restA(ti),rest_invEdgeMat[ti]);
    }

    return energyList.sum();
}

double Total_Generalized_Content_3D::compute_total_generalized_content(const Eigen::Matrix3Xd &vertices,
                                                                       Eigen::VectorXd &generalized_content_list) const {
//    generalized_content_list = Eigen::VectorXd::Zero(T.cols());
    generalized_content_list.setZero(T.cols());
    int vDim = 3;
    int simplex_size = 4; //tetrahedron

#pragma omp parallel
#pragma omp for
    for (int i = 0; i < free_tetI.size(); ++i) {
        auto ti = free_tetI(i);
        generalized_content_list(ti) = compute_generalized_TetVolume(vertices.col(T(0,ti)),
                                                                     vertices.col(T(1,ti)),
                                                                     vertices.col(T(2,ti)),
                                                                     vertices.col(T(3,ti)),
                                                                     squared_restA(ti),rest_invEdgeMat[ti]);
    }

    return generalized_content_list.sum();
}

double Total_Generalized_Content_3D::compute_total_generalized_negative_content(const Eigen::Matrix3Xd &vertices,
                                                                       Eigen::VectorXd &generalized_negative_content_list) const {
//    generalized_content_list = Eigen::VectorXd::Zero(T.cols());
    generalized_negative_content_list.setZero(T.cols());
    int vDim = 3;
    int simplex_size = 4; //tetrahedron

#pragma omp parallel
#pragma omp for
    for (int i = 0; i < free_tetI.size(); ++i) {
        auto ti = free_tetI(i);
        generalized_negative_content_list(ti) = compute_generalized_negative_TetVolume(vertices.col(T(0,ti)),
                                                                     vertices.col(T(1,ti)),
                                                                     vertices.col(T(2,ti)),
                                                                     vertices.col(T(3,ti)),
                                                                     squared_restA(ti),rest_invEdgeMat[ti]);
    }

    return generalized_negative_content_list.sum();
}

double
Total_Generalized_Content_3D::compute_generalized_TetVolume(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
                                                            const Eigen::Vector3d &v3, const Eigen::Vector3d &v4,
                                                            double squared_rest_volume,
                                                            const Eigen::Matrix3d &rest_inverse_EdgeMat) const {
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

double
Total_Generalized_Content_3D::compute_generalized_negative_TetVolume(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
                                                            const Eigen::Vector3d &v3, const Eigen::Vector3d &v4,
                                                            double squared_rest_volume,
                                                            const Eigen::Matrix3d &rest_inverse_EdgeMat) const {
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
                ) - volume;
}

double Total_Generalized_Content_3D::compute_total_generalized_content_with_gradient(const Eigen::Matrix3Xd &vertices,
                                                                                     Eigen::Matrix3Xd &grad) const {
    int vDim = 3;
    int simplex_size = 4; // tetrahedron
//    double energy = 0.0;
    Eigen::VectorXd energyList = Eigen::VectorXd::Zero(free_tetI.size());
    std::vector<Eigen::Matrix3Xd> gradList(free_tetI.size());
    grad = Eigen::Matrix3Xd::Zero(vDim, vertices.cols());

#pragma omp parallel
#pragma omp for
    for (auto i = 0; i < free_tetI.size(); ++i) {
        auto fi = free_tetI(i);
        int i1, i2, i3, i4;
        i1 = T(0, fi);
        i2 = T(1, fi);
        i3 = T(2, fi);
        i4 = T(3, fi);

//        Eigen::Matrix3Xd g;
        energyList(i) = compute_generalized_TetVolume_with_gradient(
                vertices.col(i1), vertices.col(i2), vertices.col(i3), vertices.col(i4),
                squared_restA(fi),rest_invEdgeMat[fi],
                pFpx_list[fi],gradList[i]);

//        grad.col(i1) += g.col(0);
//        grad.col(i2) += g.col(1);
//        grad.col(i3) += g.col(2);
//        grad.col(i4) += g.col(3);
    }

    // get gradient
    for (int i = 0; i < free_tetI.size(); i++)
    {
        int fi = free_tetI(i);
        grad.col(T(0, fi)) += gradList[i].col(0);
        grad.col(T(1, fi)) += gradList[i].col(1);
        grad.col(T(2, fi)) += gradList[i].col(2);
        grad.col(T(3, fi)) += gradList[i].col(3);
    }

    return energyList.sum();
}

double Total_Generalized_Content_3D::compute_generalized_TetVolume_with_gradient(const Eigen::Vector3d &v1,
                                                                                 const Eigen::Vector3d &v2,
                                                                                 const Eigen::Vector3d &v3,
                                                                                 const Eigen::Vector3d &v4,
                                                                                 double squared_rest_volume,
                                                                                 const Eigen::Matrix3d &rest_inverse_EdgeMat,
                                                                                 const Eigen::MatrixXd &pFpx,
                                                                                 Eigen::Matrix3Xd &grad) const {
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

//    grad.resize(vert.rows(), vert.cols());
    grad.resize(3, 4);
    grad = energy_constant_1 * volume * dV_dx;
    Eigen::VectorXd px = pFpx.transpose() * (grad_constant_1 * squared_rest_volume * vec_pF);
    grad.col(0) += px.segment<3>(0);
    grad.col(1) += px.segment<3>(3);
    grad.col(2) += px.segment<3>(6);
    grad.col(3) += px.segment<3>(9);
    grad /= energy;

    return energy;
}

double Total_Generalized_Content_3D::compute_total_generalized_content_with_gradient_and_sTGC_projectedHessian(
        const Eigen::Matrix3Xd &vertices, const Eigen::VectorXi &freeI, const Eigen::Matrix4Xi &T_free,
        Eigen::VectorXd &generalized_content_list, Eigen::Matrix3Xd &grad, SpMat &Hess) const {
    int vDim = 3;
    int simplex_size = 4;
    generalized_content_list.resize(T.cols());
    std::vector<Eigen::Matrix3Xd> gradList(T.cols());
    grad = Eigen::Matrix3Xd::Zero(3, vertices.cols());

    std::vector<Eigen::Triplet<double>> tripletList(simplex_size * simplex_size * vDim * vDim * T.cols());

#pragma omp parallel
#pragma omp for
    for (auto i = 0; i < T.cols(); ++i) {
        Eigen::MatrixXd  hess;
        generalized_content_list(i) = compute_generalized_TetVolume_with_gradient_projected_subtracted_Hessian(
                vertices.col(T(0,i)), vertices.col(T(1,i)),
                vertices.col(T(2,i)), vertices.col(T(3,i)),
                squared_restA(i), restA(i),
                rest_invEdgeMat[i],
                pFpx_list[i],
                gradList[i], hess);

        // update Hessian of free vertices
        int current_index = i * simplex_size * simplex_size * vDim * vDim;
        Eigen::Vector4i indices = T_free.col(i);
        for (int j = 0; j < simplex_size; ++j) {
            int idx_j = indices(j);
            for (int k = 0; k < simplex_size; ++k) {
                int idx_k = indices(k);
                if (idx_j != -1 && idx_k != -1) {
                    for (int l = 0; l < vDim; ++l) {
                        for (int n = 0; n < vDim; ++n) {
                            tripletList[current_index] = Eigen::Triplet<double>(idx_j * vDim + l, idx_k * vDim + n,
                                                                                hess(j * vDim + l, k * vDim + n));
                            ++current_index;
                        }
                    }
                }
            }
        }

    }

    // get gradient
    for (int i = 0; i < T.cols(); i++)
    {
        grad.col(T(0, i)) += gradList[i].col(0);
        grad.col(T(1, i)) += gradList[i].col(1);
        grad.col(T(2, i)) += gradList[i].col(2);
        grad.col(T(3, i)) += gradList[i].col(3);
    }

    // add small positive values to the diagonal of Hessian
    for (auto i = 0; i < vDim * freeI.size(); ++i) {
        tripletList.emplace_back(i, i, 1e-8);
    }

    // get Hessian on free vertices
    Hess.resize(vDim * freeI.size(), vDim * freeI.size());
    Hess.setFromTriplets(tripletList.begin(), tripletList.end());

    return generalized_content_list.sum();
}

double Total_Generalized_Content_3D::compute_generalized_TetVolume_with_gradient_projected_subtracted_Hessian(
        const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d &v3, const Eigen::Vector3d &v4,
        double squared_rest_volume, double rest_volume,
        const Eigen::Matrix3d &rest_inverse_EdgeMat, const Eigen::MatrixXd &pFpx,
        Eigen::Matrix3Xd &grad, Eigen::MatrixXd &Hess) const {
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
    Vector9d vec_pF;
    vectorize3x3(pF, vec_pF);

//    grad.resize(vert.rows(), vert.cols());
    grad.resize(3, 4);
    grad = energy_constant_1 * volume * dV_dx;
    Eigen::VectorXd px = pFpx.transpose() * (grad_constant_1 * squared_rest_volume * vec_pF);
    grad.col(0) += px.segment<3>(0);
    grad.col(1) += px.segment<3>(3);
    grad.col(2) += px.segment<3>(6);
    grad.col(3) += px.segment<3>(9);
    grad /= energy;

    // PSD projected Hessian of (generalized content - signed content)
    Eigen::Matrix3d U, V;
    Eigen::Vector3d singular_values;
    compute_SVD(deformation_gradient, U, singular_values, V);
    double s1 = singular_values[0];
    double s2 = singular_values[1];
    double s3 = singular_values[2];
    // S-centric invariants
//    double I1 = s1 + s2 + s3;
    double I2 = s1*s1 + s2*s2 + s3*s3;
    double I3 = s1*s2*s3;
    double s11 = s1*s1;
    double s22 = s2*s2;
    double s33 = s3*s3;
    double s12 = s1*s2;
    double s13 = s1*s3;
    double s23 = s2*s3;
    double II_C_star = s12 * s12 + s13 * s13 + s23 * s23;
    // Analytic Eigensystem of f-Hessian
    Eigen::Matrix<double,9,9> f_Hess;
    f_Hess.setZero();
    double energy_density = energy / rest_volume;
    Vector9d eigen_vec;
    // twist 1
    double eigen_value_twist_1 = (energy_constant_2 + 2 * energy_constant_3 * (s11 + s23) + energy_constant_1 * s1 * I3)/energy_density - s1;
    if (eigen_value_twist_1 > 0) {
//        Vector9d eigen_vec_twist_1;
        // 0 0 0
        // 0 0 1
        // 0 -1 0
        vectorize3x3((U.col(1) * V.col(2).transpose() - U.col(2) * V.col(1).transpose()) / sqrt(2.), eigen_vec);
        f_Hess += eigen_value_twist_1 * eigen_vec * eigen_vec.transpose();
    }
    // twist 2
    double eigen_value_twist_2 = (energy_constant_2 + 2 * energy_constant_3 * (s22 + s13) + energy_constant_1 * s2 * I3)/energy_density - s2;
    if (eigen_value_twist_2 > 0) {
//        Vector9d eigen_vec_twist_2;
        // 0 0 1
        // 0 0 0
        // -1 0 0
        vectorize3x3((U.col(0) * V.col(2).transpose() - U.col(2) * V.col(0).transpose()) / sqrt(2.), eigen_vec);
        f_Hess += eigen_value_twist_2 * eigen_vec * eigen_vec.transpose();
    }
    // twist 3
    double eigen_value_twist_3 = (energy_constant_2 + 2 * energy_constant_3 * (s33 + s12) + energy_constant_1 * s3 * I3)/energy_density - s3;
    if (eigen_value_twist_3 > 0) {
//        Vector9d eigen_vec_twist_3;
        // 0 -1 0
        // 1 0 0
        // 0 0 0
        vectorize3x3((U.col(1) * V.col(0).transpose() - U.col(0) * V.col(1).transpose()) / sqrt(2.), eigen_vec);
        f_Hess += eigen_value_twist_3 * eigen_vec * eigen_vec.transpose();
    }
    // flip 1
    double eigen_value_flip_1 = (energy_constant_2 + 2 * energy_constant_3 * (s11 - s23) - energy_constant_1 * s1 * I3)/energy_density + s1;
    if (eigen_value_flip_1 > 0) {
//        Vector9d eigen_vec_flip_1;
        // 0 0 0
        // 0 0 1
        // 0 1 0
        vectorize3x3((U.col(1) * V.col(2).transpose() + U.col(2) * V.col(1).transpose()) / sqrt(2.), eigen_vec);
        f_Hess += eigen_value_flip_1 * eigen_vec * eigen_vec.transpose();
    }
    // flip 2
    double eigen_value_flip_2 = (energy_constant_2 + 2 * energy_constant_3 * (s22 - s13) - energy_constant_1 * s2 * I3)/energy_density + s2;
    if (eigen_value_flip_2 > 0) {
//        Vector9d eigen_vec_flip_2;
        // 0 0 1
        // 0 0 0
        // 1 0 0
        vectorize3x3((U.col(0) * V.col(2).transpose() + U.col(2) * V.col(0).transpose()) / sqrt(2.), eigen_vec);
        f_Hess += eigen_value_flip_2 * eigen_vec * eigen_vec.transpose();
    }
    // flip 3
    double eigen_value_flip_3 = (energy_constant_2 + 2 * energy_constant_3 * (s33 - s12) - energy_constant_1 * s3 * I3)/energy_density + s3;
    if (eigen_value_flip_3 > 0) {
//        Vector9d eigen_vec_flip_3;
        // 0 1 0
        // 1 0 0
        // 0 0 0
        vectorize3x3((U.col(0) * V.col(1).transpose() + U.col(1) * V.col(0).transpose()) / sqrt(2.), eigen_vec);
        f_Hess += eigen_value_flip_3 * eigen_vec * eigen_vec.transpose();
    }
    //// scale
    double energy_cube = energy_density * energy_density * energy_density;
    double a11 = (energy_constant_4 + 2 * energy_constant_3 * s22 * s33 + energy_constant_2 * (s22 + s33))
            * (energy_constant_2 + energy_constant_1 * s22 * s33 + 2 * energy_constant_3 * (s22 + s33)) / energy_cube;
    double a22 = (energy_constant_4 + 2 * energy_constant_3 * s11 * s33 + energy_constant_2 * (s11 + s33))
                 * (energy_constant_2 + energy_constant_1 * s11 * s33 + 2 * energy_constant_3 * (s11 + s33)) / energy_cube;
    double a33 = (energy_constant_4 + 2 * energy_constant_3 * s11 * s22 + energy_constant_2 * (s11 + s22))
                 * (energy_constant_2 + energy_constant_1 * s11 * s22 + 2 * energy_constant_3 * (s11 + s22)) / energy_cube;

    double a12 = (energy_constant_2 * (energy_constant_1 * I3 * s3 * (I2 + s33) + 2 * energy_constant_3 * s12 * (I2 - s33) - energy_density * I2 * s3)
                  + 4 * energy_constant_3 * energy_constant_3 * s12 * (II_C_star - s33 * s33)
                  + 2 * energy_constant_3 * (2 * energy_constant_4 * s12 - energy_density * II_C_star * s3 + energy_constant_1 * I3 * (II_C_star + s12 * s12) * s3)
                  + (energy_constant_1 * I3 * I3 * (energy_constant_1 * I3 - energy_density) - energy_constant_4 * (energy_density - 2 * energy_constant_1 * I3)) * s3
                  - energy_constant_2 * energy_constant_2 * s12 * s12 ) / energy_cube;
    double a13 = (energy_constant_2 * (energy_constant_1 * I3 * s2 * (I2 + s22) + 2 * energy_constant_3 * s13 * (I2 - s22) - energy_density * I2 * s2)
            + 4 * energy_constant_3 * energy_constant_3 * s13 * (II_C_star - s22 * s22)
            + 2 * energy_constant_3 * (2 * energy_constant_4 * s13 - energy_density * II_C_star * s2 + energy_constant_1 * I3 * (II_C_star + s13 * s13) * s2)
            + (energy_constant_1 * I3 * I3 * (energy_constant_1 * I3 - energy_density) - energy_constant_4 * (energy_density - 2 * energy_constant_1 * I3)) * s2
            - energy_constant_2 * energy_constant_2 * s13 * s13 ) / energy_cube;
    double a23 = (energy_constant_2 * (energy_constant_1 * I3 * s1 * (I2 + s11) + 2 * energy_constant_3 * s23 * (I2 - s11) - energy_density * I2 * s1)
                  + 4 * energy_constant_3 * energy_constant_3 * s23 * (II_C_star - s11 * s11)
                  + 2 * energy_constant_3 * (2 * energy_constant_4 * s23 - energy_density * II_C_star * s1 + energy_constant_1 * I3 * (II_C_star + s23 * s23) * s1)
                  + (energy_constant_1 * I3 * I3 * (energy_constant_1 * I3 - energy_density) - energy_constant_4 * (energy_density - 2 * energy_constant_1 * I3)) * s1
                  - energy_constant_2 * energy_constant_2 * s23 * s23 ) / energy_cube;

    if (a12 == 0 && a13 == 0 && a23 == 0) {
//        eigen_value_scale_1 = a11;
//        eigen_value_scale_2 = a22;
//        eigen_value_scale_3 = a33;
//        eigen_vec_scale_1 = vec_d1;
//        eigen_vec_scale_2 = vec_d2;
//        eigen_vec_scale_3 = vec_d3;
        if (a11 > 0) {
            vectorize3x3(U.col(0) * V.col(0).transpose(), eigen_vec);
            f_Hess += a11 * eigen_vec * eigen_vec.transpose();
        }
        if (a22 > 0) {
            vectorize3x3(U.col(1) * V.col(1).transpose(), eigen_vec);
            f_Hess += a22 * eigen_vec * eigen_vec.transpose();
        }
        if (a33 > 0) {
            vectorize3x3(U.col(2) * V.col(2).transpose(), eigen_vec);
            f_Hess += a33 * eigen_vec * eigen_vec.transpose();
        }
    }
    else {
        Vector9d vec_d1, vec_d2, vec_d3;
        vectorize3x3(U.col(0) * V.col(0).transpose(), vec_d1);
        vectorize3x3(U.col(1) * V.col(1).transpose(), vec_d2);
        vectorize3x3(U.col(2) * V.col(2).transpose(), vec_d3);
        Eigen::Matrix3d matA;
        matA << a11, a12, a13,
                a12, a22, a23,
                a13, a23, a33;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(matA);
        Eigen::Vector3d matA_eigenVals = eigenSolver.eigenvalues();
        // eigen vectors has been normalized to have length 1
        Eigen::Matrix3d matA_eigenVecs = eigenSolver.eigenvectors();
        //
        if (matA_eigenVals[0] > 0) {
            eigen_vec = matA_eigenVecs(0,0) * vec_d1 + matA_eigenVecs(1,0) * vec_d2 + matA_eigenVecs(2,0) * vec_d3;
            f_Hess += matA_eigenVals[0] * eigen_vec * eigen_vec.transpose();
        }
        if (matA_eigenVals[1] > 0) {
            eigen_vec = matA_eigenVecs(0,1) * vec_d1 + matA_eigenVecs(1,1) * vec_d2 + matA_eigenVecs(2,1) * vec_d3;
            f_Hess += matA_eigenVals[1] * eigen_vec * eigen_vec.transpose();
        }
        if (matA_eigenVals[2] > 0) {
            eigen_vec = matA_eigenVecs(0,2) * vec_d1 + matA_eigenVecs(1,2) * vec_d2 + matA_eigenVecs(2,2) * vec_d3;
            f_Hess += matA_eigenVals[2] * eigen_vec * eigen_vec.transpose();
        }
    }
    // PSD x-Hessian
    Hess = pFpx.transpose() * (rest_volume * f_Hess) * pFpx;

    return energy;
}
