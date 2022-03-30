//
// Created by Charles Du on 3/28/22.
//

#include "Total_Generalized_Content.h"

#include <Eigen/Eigenvalues>
#include "deformation_gradient_util.h"


void Total_Generalized_Content::initialize(const Eigen::MatrixXd &rest_vertices, const Eigen::Matrix3Xi &faces,
                                           const std::string &form, double alpha, double lambda1, double lambda2,
                                           double k) {
    //
    F = faces;
    param_alpha = alpha;
    param_lambda1 = lambda1;
    param_lambda2 = lambda2;
    param_k = k;
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
        // compute inverse edge matrix and pFpx
        Eigen::Vector2d v1(0,0);
        double regular_tri_edge_length = sqrt(regular_tri_squared_length);
        Eigen::Vector2d v2(regular_tri_edge_length, 0);
        Eigen::Vector2d v3(regular_tri_edge_length/2, sqrt(3.)*regular_tri_edge_length/2);
        Eigen::Matrix2d edgeMat;
        compute_edge_matrix(v1,v2,v3, edgeMat);
        Eigen::Matrix2d invEdgeMat = edgeMat.inverse();
        Eigen::MatrixXd pFpx;
        compute_flatten_pFpx(invEdgeMat, pFpx);
        rest_invEdgeMat.resize(F.cols(), invEdgeMat);
        pFpx_list.resize(F.cols(), pFpx);
    } else { // harmonic form
        rest_invEdgeMat.reserve(F.cols());
        Eigen::Matrix2d rest_edge_mat;
        for (int i = 0; i < F.cols(); ++i) {
            rest_invEdgeMat.emplace_back();
            compute_edge_matrix(rest_vertices.col(F(0,i)),
                                rest_vertices.col(F(1,i)),
                                rest_vertices.col(F(2,i)),
                                rest_edge_mat);
            rest_invEdgeMat.back() = rest_edge_mat.inverse();
        }
        // compute derivative of deformation gradient w.r.t. target vertices
        pFpx_list.reserve(F.cols());
        for (int i = 0; i < F.cols(); ++i) {
            pFpx_list.emplace_back();
            compute_flatten_pFpx(rest_invEdgeMat[i], pFpx_list.back());
        }
    }

    // compute scaled squared areas of auxiliary triangles: (2*alpha*lambda2*k^2 + alpha^2)*(A_rest)^2
    scaled_squared_restA.resize(F.cols());
    for (int i = 0; i < F.cols(); ++i) {
        scaled_squared_restA(i) = restA(i) * restA(i);
    }
    scaled_squared_restA *= (2*alpha*lambda2*k*k + alpha*alpha);

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
    scaled_Dirichlet_coef *= (alpha * lambda1/4);

    squared_targetA_scale_coeff = 1 + 2*alpha*lambda2;

    // coefficients for analytic eigen system of Hessian
    two_alpha_lambda1 = 2 * alpha * lambda1;
    one_plus_two_alpha_lambda2 = 1. + 2 * alpha * lambda2;
    coeff_diag = alpha * (alpha + 2 * lambda2 * k * k);
    coeff_off_diag = alpha*alpha*(2 - 4*lambda1*lambda1 + 4*(1+k*k)*alpha*lambda2 + 8*k*k*lambda2*lambda2);
    coeff_off_diag_subtracted = 2*alpha*(alpha - 2*alpha*lambda1*lambda1 + 2*k*k*lambda2 + 2*alpha*alpha*lambda2 + 4*k*k*alpha*lambda2*lambda2);
}

double Total_Generalized_Content::compute_total_generalized_content(const Eigen::Matrix2Xd &vertices) const {
    Eigen::VectorXd energyList = Eigen::VectorXd::Zero(F.cols());
    int vDim = 2;
    int simplex_size = 3; //triangle

    for (int i = 0; i < free_faceI.size(); ++i) {
        auto fi = free_faceI(i);
        Eigen::Matrix2Xd vert(vDim, simplex_size);
        vert.col(0) = vertices.col(F(0, fi));
        vert.col(1) = vertices.col(F(1, fi));
        vert.col(2) = vertices.col(F(2, fi));

        energyList(fi) = compute_generalized_TriArea(vert, scaled_Dirichlet_coef.col(fi),scaled_squared_restA(fi));
    }

    return energyList.sum();
}

double Total_Generalized_Content::compute_generalized_TriArea(const Eigen::Matrix2Xd &vert,
                                                              const Eigen::Vector3d &Dirichlet_coef,
                                                              double scaled_squared_rest_area) const {
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto e1 = v2 - v3;
    auto e2 = v3 - v1;
    auto e3 = v1 - v2;

    double area = compute_tri_signed_area(v1,v2,v3);

    return sqrt(squared_targetA_scale_coeff*area*area + Dirichlet_coef(0) * e1.squaredNorm()
                + Dirichlet_coef(1) * e2.squaredNorm() + Dirichlet_coef(2) * e3.squaredNorm()
                + scaled_squared_rest_area
    );
}

double Total_Generalized_Content::compute_total_generalized_content(const Eigen::Matrix2Xd &vertices,
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

        energyList(fi) = compute_generalized_TriArea(
                vert, scaled_Dirichlet_coef.col(fi),scaled_squared_restA(fi));
    }

    return energyList.sum();
}

double Total_Generalized_Content::compute_total_generalized_content_with_gradient(const Eigen::Matrix2Xd &vertices,
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
        energy += compute_generalized_TriArea_with_gradient(
                vert,scaled_Dirichlet_coef.col(fi),scaled_squared_restA(fi),g);

        grad.col(i1) += g.col(0);
        grad.col(i2) += g.col(1);
        grad.col(i3) += g.col(2);
    }

    return energy;
}

double Total_Generalized_Content::compute_generalized_TriArea_with_gradient(const Eigen::Matrix2Xd &vert,
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
    double energy = sqrt(squared_targetA_scale_coeff*area*area + d1 * e1.squaredNorm()
                         + d2 * e2.squaredNorm() + d3 * e3.squaredNorm()
                         + scaled_squared_rest_area);

    // partial(2*area)/partial(v1,v2,v3)
    auto h1 = rotate_90deg(e1);
    auto h2 = rotate_90deg(e2);
    auto h3 = rotate_90deg(e3);

    //
    double f1 = (0.5+param_alpha*param_lambda2) * area;

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

double Total_Generalized_Content::compute_total_generalized_content_with_gradient_and_sTGC_projectedHessian(
        const Eigen::Matrix2Xd &vertices, const Eigen::VectorXi &freeI, const Eigen::Matrix3Xi &F_free,
        Eigen::VectorXd &generalized_content_list, Eigen::Matrix2Xd &grad, SpMat &Hess) const {
    int vDim = 2;
    generalized_content_list.resize(F.cols());
    grad = Eigen::Matrix2Xd::Zero(2, vertices.cols());

    std::vector<Eigen::Triplet<double>> tripletList(3 * 3 * vDim * vDim * F.cols());

#pragma omp parallel
#pragma omp for
    for (auto i = 0; i < F.cols(); ++i) {
        int i1, i2, i3;
        i1 = F(0, i);
        i2 = F(1, i);
        i3 = F(2, i);

        Eigen::MatrixXd vert(vDim, 3);
        vert.col(0) = vertices.col(i1);
        vert.col(1) = vertices.col(i2);
        vert.col(2) = vertices.col(i3);
        Eigen::Vector3d r = restD.col(i);


        Eigen::Matrix2Xd g;
        Eigen::MatrixXd  hess;
        generalized_content_list(i) = compute_generalized_TriArea_with_gradient_projected_subtracted_Hessian(vert,
                                                                                                 scaled_Dirichlet_coef.col(i),
                                                                                                 scaled_squared_restA(i),
                                                                                                 restA(i),
                                                                                                 rest_invEdgeMat[i],
                                                                                                 pFpx_list[i],
                                                                                                 two_alpha_lambda1,
                                                                                                 one_plus_two_alpha_lambda2,
                                                                                                 coeff_diag,
                                                                                                 coeff_off_diag_subtracted,
                                                                                                 g, hess);

#pragma omp critical
        {
            grad.col(i1) += g.col(0);
            grad.col(i2) += g.col(1);
            grad.col(i3) += g.col(2);
        }

        // update Hessian of free vertices
        int current_index = i * 3 * 3 * vDim * vDim;
        Eigen::Vector3i indices = F_free.col(i);
        for (int j = 0; j < 3; ++j) {
            int idx_j = indices(j);
            for (int k = 0; k < 3; ++k) {
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

    // add small positive values to the diagonal of Hessian
    for (auto i = 0; i < vDim * freeI.size(); ++i) {
        tripletList.emplace_back(i, i, 1e-8);
    }

    // get Hessian on free vertices
    Hess.resize(vDim * freeI.size(), vDim * freeI.size());
    Hess.setFromTriplets(tripletList.begin(), tripletList.end());

    return generalized_content_list.sum();
}

double
Total_Generalized_Content::compute_generalized_TriArea_with_gradient_projected_subtracted_Hessian(const Eigen::Matrix2Xd &vert,
                                                                                      const Eigen::Vector3d &Dirichlet_coef,
                                                                                      double scaled_squared_rest_area,
                                                                                      double rest_area,
                                                                                      const Eigen::Matrix2d &rest_inverse_EdgeMat,
                                                                                      const Eigen::MatrixXd &pFpx,
                                                                                      double two_alpha_lambda1,
                                                                                      double one_plus_two_alpha_lambda2,
                                                                                      double coeff_diag,
                                                                                      double coeff_off_diag_subtracted,
                                                                                      Eigen::Matrix2Xd &grad,
                                                                                      Eigen::MatrixXd &Hess) const {
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
    double energy = sqrt(squared_targetA_scale_coeff*area*area + d1 * e1.squaredNorm()
                         + d2 * e2.squaredNorm() + d3 * e3.squaredNorm()
                         + scaled_squared_rest_area);

    // partial(2*area)/partial(v1,v2,v3)
    auto h1 = rotate_90deg(e1);
    auto h2 = rotate_90deg(e2);
    auto h3 = rotate_90deg(e3);

    //
    double f1 = (0.5+param_alpha*param_lambda2) * area;

    auto d1e1 = d1 * e1;
    auto d2e2 = d2 * e2;
    auto d3e3 = d3 * e3;

    grad.resize(vert.rows(), vert.cols());
    grad.col(0) = f1 * h1 + (d2e2 - d3e3);
    grad.col(1) = f1 * h2 + (d3e3 - d1e1);
    grad.col(2) = f1 * h3 + (d1e1 - d2e2);

    grad /= energy;

    // PSD projected Hessian of (generalized triangle area - signed triangle area)
    //
    Eigen::Matrix2d target_edge_matrix;
    compute_edge_matrix(v1, v2, v3, target_edge_matrix);
    Eigen::Matrix2d deformation_gradient = target_edge_matrix * rest_inverse_EdgeMat;
    Eigen::Matrix2d U, V;
    Eigen::Vector2d singular_values;
    compute_SVD(deformation_gradient, U, singular_values, V);
    double s1 = singular_values[0];
    double s2 = singular_values[1];
    // S-centric invariants
    double I1 = s1 + s2;
    double I2 = s1*s1 + s2*s2;
    double I3 = s1*s2;
    // Analytic Eigensystem of f-Hessian
    double energy_density = energy / rest_area;
    //// twist
    double eigen_value_twist = (two_alpha_lambda1 + I3 * one_plus_two_alpha_lambda2)/energy_density - 1;
    Eigen::VectorXd eigen_vec_twist;
    vectorize((U.col(1) * V.col(0).transpose() - U.col(0) * V.col(1).transpose())/sqrt(2.),
              eigen_vec_twist);
    //// flip
    double eigen_value_flip = (two_alpha_lambda1 - I3 * one_plus_two_alpha_lambda2)/energy_density + 1;
    Eigen::VectorXd eigen_vec_flip;
    vectorize((U.col(1) * V.col(0).transpose() + U.col(0) * V.col(1).transpose())/sqrt(2.),
              eigen_vec_flip);
    //// scale
    Eigen::VectorXd vec_d1, vec_d2;
    vectorize(U.col(0) * V.col(0).transpose(), vec_d1);
    vectorize(U.col(1) * V.col(1).transpose(), vec_d2);

    double energy_cube = energy_density * energy_density * energy_density;
    double a11 = (coeff_diag + two_alpha_lambda1*s2*s2)*(two_alpha_lambda1+one_plus_two_alpha_lambda2*s2*s2)/energy_cube;
    double a22 = (coeff_diag + two_alpha_lambda1*s1*s1)*(two_alpha_lambda1+one_plus_two_alpha_lambda2*s1*s1)/energy_cube;
    double a12 = (coeff_off_diag_subtracted*I3 + one_plus_two_alpha_lambda2*I3*(one_plus_two_alpha_lambda2*I3*I3+two_alpha_lambda1*I2))/energy_cube-1;
    Eigen::Matrix2d matA;
    matA << a11, a12,
            a12, a22;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenSolver(matA);
    auto matA_eigenVals = eigenSolver.eigenvalues();
    double eigen_value_scale1 = matA_eigenVals[0];
    double eigen_value_scale2 = matA_eigenVals[1];
    double beta = (eigen_value_scale1 - a22)/a12;
    double norm_beta_1 = sqrt(1. + beta*beta);
    Eigen::VectorXd eigen_vec_scale1 = (beta * vec_d1 + vec_d2) / norm_beta_1;
    Eigen::VectorXd eigen_vec_scale2 = (vec_d1 - beta * vec_d2) / norm_beta_1;
    // PSD f-Hessian
    Eigen::Matrix4d f_Hess;
    f_Hess.setZero();
    if (eigen_value_twist > 0) {
        f_Hess += eigen_value_twist * eigen_vec_twist * eigen_vec_twist.transpose();
    }
    if (eigen_value_flip > 0) {
        f_Hess += eigen_value_flip * eigen_vec_flip * eigen_vec_flip.transpose();
    }
    if (eigen_value_scale1 > 0) {
        f_Hess += eigen_value_scale1 * eigen_vec_scale1 * eigen_vec_scale1.transpose();
    }
    if (eigen_value_scale2 > 0) {
        f_Hess += eigen_value_scale2 * eigen_vec_scale2 * eigen_vec_scale2.transpose();
    }
    // PSD x-Hessian
    Hess = pFpx.transpose() * (rest_area * f_Hess) * pFpx;

    return energy;
}

double Total_Generalized_Content::compute_total_generalized_content_with_gradient_and_TGC_projectedHessian(
        const Eigen::Matrix2Xd &vertices, const Eigen::VectorXi &freeI, const Eigen::Matrix3Xi &F_free,
        Eigen::VectorXd &generalized_content_list, Eigen::Matrix2Xd &grad, SpMat &Hess) const {
    int vDim = 2;
    generalized_content_list.resize(F.cols());
    grad = Eigen::Matrix2Xd::Zero(2, vertices.cols());

    std::vector<Eigen::Triplet<double>> tripletList(3 * 3 * vDim * vDim * F.cols());

#pragma omp parallel
#pragma omp for
    for (auto i = 0; i < F.cols(); ++i) {
        int i1, i2, i3;
        i1 = F(0, i);
        i2 = F(1, i);
        i3 = F(2, i);

        Eigen::MatrixXd vert(vDim, 3);
        vert.col(0) = vertices.col(i1);
        vert.col(1) = vertices.col(i2);
        vert.col(2) = vertices.col(i3);
        Eigen::Vector3d r = restD.col(i);


        Eigen::Matrix2Xd g;
        Eigen::MatrixXd  hess;
        generalized_content_list(i) = compute_generalized_TriArea_with_gradient_projectedHessian(vert,
                                                                                                 scaled_Dirichlet_coef.col(i),
                                                                                                 scaled_squared_restA(i),
                                                                                                 restA(i),
                                                                                                 rest_invEdgeMat[i],
                                                                                                 pFpx_list[i],
                                                                                                 two_alpha_lambda1,
                                                                                                 one_plus_two_alpha_lambda2,
                                                                                                 coeff_diag,
                                                                                                 coeff_off_diag,
                                                                                                 g, hess);

#pragma omp critical
        {
            grad.col(i1) += g.col(0);
            grad.col(i2) += g.col(1);
            grad.col(i3) += g.col(2);
        }

        // update Hessian of free vertices
        int current_index = i * 3 * 3 * vDim * vDim;
        Eigen::Vector3i indices = F_free.col(i);
        for (int j = 0; j < 3; ++j) {
            int idx_j = indices(j);
            for (int k = 0; k < 3; ++k) {
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

    // add small positive values to the diagonal of Hessian
    for (auto i = 0; i < vDim * freeI.size(); ++i) {
        tripletList.emplace_back(i, i, 1e-8);
    }

    // get Hessian on free vertices
    Hess.resize(vDim * freeI.size(), vDim * freeI.size());
    Hess.setFromTriplets(tripletList.begin(), tripletList.end());

    return generalized_content_list.sum();
}

double
Total_Generalized_Content::compute_generalized_TriArea_with_gradient_projectedHessian(const Eigen::Matrix2Xd &vert,
                                                                                      const Eigen::Vector3d &Dirichlet_coef,
                                                                                      double scaled_squared_rest_area,
                                                                                      double rest_area,
                                                                                      const Eigen::Matrix2d &rest_inverse_EdgeMat,
                                                                                      const Eigen::MatrixXd &pFpx,
                                                                                      double two_alpha_lambda1,
                                                                                      double one_plus_two_alpha_lambda2,
                                                                                      double coeff_diag,
                                                                                      double coeff_off_diag,
                                                                                      Eigen::Matrix2Xd &grad,
                                                                                      Eigen::MatrixXd &Hess) const {
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
    double energy = sqrt(squared_targetA_scale_coeff*area*area + d1 * e1.squaredNorm()
                         + d2 * e2.squaredNorm() + d3 * e3.squaredNorm()
                         + scaled_squared_rest_area);

    // partial(2*area)/partial(v1,v2,v3)
    auto h1 = rotate_90deg(e1);
    auto h2 = rotate_90deg(e2);
    auto h3 = rotate_90deg(e3);

    //
    double f1 = (0.5+param_alpha*param_lambda2) * area;

    auto d1e1 = d1 * e1;
    auto d2e2 = d2 * e2;
    auto d3e3 = d3 * e3;

    grad.resize(vert.rows(), vert.cols());
    grad.col(0) = f1 * h1 + (d2e2 - d3e3);
    grad.col(1) = f1 * h2 + (d3e3 - d1e1);
    grad.col(2) = f1 * h3 + (d1e1 - d2e2);

    grad /= energy;

    // PSD projected Hessian of generalized triangle area
    // deformation gradient and singular values
    Eigen::Matrix2d target_edge_matrix;
    compute_edge_matrix(v1, v2, v3, target_edge_matrix);
    Eigen::Matrix2d deformation_gradient = target_edge_matrix * rest_inverse_EdgeMat;
    Eigen::Matrix2d U, V;
    Eigen::Vector2d singular_values;
    compute_SVD(deformation_gradient, U, singular_values, V);
    double s1 = singular_values[0];
    double s2 = singular_values[1];
    // S-centric invariants
    double I1 = s1 + s2;
    double I2 = s1*s1 + s2*s2;
    double I3 = s1*s2;
    // Analytic Eigensystem of f-Hessian
    double energy_density = energy / rest_area;
    //// twist
    double eigen_value_twist = (two_alpha_lambda1 + I3 * one_plus_two_alpha_lambda2)/energy_density;
    Eigen::VectorXd eigen_vec_twist;
    vectorize((U.col(1) * V.col(0).transpose() - U.col(0) * V.col(1).transpose())/sqrt(2.),
              eigen_vec_twist);
    //// flip
    double eigen_value_flip = (two_alpha_lambda1 - I3 * one_plus_two_alpha_lambda2)/energy_density;
    Eigen::VectorXd eigen_vec_flip;
    vectorize((U.col(1) * V.col(0).transpose() + U.col(0) * V.col(1).transpose())/sqrt(2.),
              eigen_vec_flip);
    //// scale
    Eigen::VectorXd vec_d1, vec_d2;
    vectorize(U.col(0) * V.col(0).transpose(), vec_d1);
    vectorize(U.col(1) * V.col(1).transpose(), vec_d2);
    double energy_cube = energy_density * energy_density * energy_density;
    //
    double a11 = (coeff_diag + two_alpha_lambda1*s2*s2)*(two_alpha_lambda1+one_plus_two_alpha_lambda2*s2*s2)/energy_cube;
    double a22 = (coeff_diag + two_alpha_lambda1*s1*s1)*(two_alpha_lambda1+one_plus_two_alpha_lambda2*s1*s1)/energy_cube;
    double a12 = (coeff_off_diag*I3 + one_plus_two_alpha_lambda2*I3*(one_plus_two_alpha_lambda2*I3*I3+two_alpha_lambda1*I2))/energy_cube;
    //
    Eigen::Matrix2d matA;
    matA << a11, a12,
            a12, a22;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenSolver(matA);
    auto matA_eigenVals = eigenSolver.eigenvalues();
    double eigen_value_scale1 = matA_eigenVals[0];
    double eigen_value_scale2 = matA_eigenVals[1];
    double beta = (eigen_value_scale1 - a22)/a12;
    double norm_beta_1 = sqrt(1. + beta*beta);
    Eigen::VectorXd eigen_vec_scale1 = (beta * vec_d1 + vec_d2) / norm_beta_1;
    Eigen::VectorXd eigen_vec_scale2 = (vec_d1 - beta * vec_d2) / norm_beta_1;
    // PSD f-Hessian
    Eigen::Matrix4d f_Hess;
    f_Hess.setZero();
    if (eigen_value_twist > 0) {
        f_Hess += eigen_value_twist * eigen_vec_twist * eigen_vec_twist.transpose();
    }
    if (eigen_value_flip > 0) {
        f_Hess += eigen_value_flip * eigen_vec_flip * eigen_vec_flip.transpose();
    }
    if (eigen_value_scale1 > 0) {
        f_Hess += eigen_value_scale1 * eigen_vec_scale1 * eigen_vec_scale1.transpose();
    }
    if (eigen_value_scale2 > 0) {
        f_Hess += eigen_value_scale2 * eigen_vec_scale2 * eigen_vec_scale2.transpose();
    }
    // PSD x-Hessian
    Hess = pFpx.transpose() * (rest_area * f_Hess) * pFpx;

    return energy;
}

