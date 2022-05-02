//
// Created by Charles Du on 4/5/22.
//

#ifndef TLC_TOTAL_GENERALIZED_CONTENT_3D_H
#define TLC_TOTAL_GENERALIZED_CONTENT_3D_H

#include "geo_util.h"

#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> SpMat;

// total smooth generalized content in 3D
class Total_Generalized_Content_3D {
public:
    Total_Generalized_Content_3D() = default;

    // initialize from rest mesh and parameters: get ready for computing TGC energy and derivatives
    // need to satisfy 6 * k^(1/3) * lambda1 + 2 * k * lambda2 = 1
    void initialize(const Eigen::Matrix3Xd &rest_vertices, const Eigen::Matrix4Xi &tets,
                    const std::string &form, double alpha, double lambda1, double lambda2, double k);

    ~Total_Generalized_Content_3D() = default;

    // compute total generalized content
    double compute_total_generalized_content(const Eigen::Matrix3Xd &vertices) const;

    // compute total generalized negative content
    double compute_total_generalized_negative_content(const Eigen::Matrix3Xd &vertices) const;

    // compute total generalized content, record generalized content of each tetrahedron
    double compute_total_generalized_content(const Eigen::Matrix3Xd &vertices,
                                             Eigen::VectorXd &generalized_content_list) const;

    // compute total generalized negative content, record generalized negative content of each tetrahedron
    double compute_total_generalized_negative_content(const Eigen::Matrix3Xd &vertices,
                                             Eigen::VectorXd &generalized_negative_content_list) const;

    // compute total generalized content and its gradient
    double compute_total_generalized_content_with_gradient(const Eigen::Matrix3Xd &vertices,
                                                           Eigen::Matrix3Xd &grad) const;


    // compute total generalized content and its gradient,
    // and PSD projected Hessian of (TGC - total signed content) on free vertices
    double compute_total_generalized_content_with_gradient_and_sTGC_projectedHessian(
            const Eigen::Matrix3Xd &vertices,
            const Eigen::VectorXi &freeI,
            const Eigen::Matrix4Xi  &T_free,
            // output
            Eigen::VectorXd &generalized_content_list,
            Eigen::Matrix3Xd &grad,
            SpMat &Hess) const;

    // compute total generalized content and its gradient,
    // and PSD projected Hessian of TGC on free vertices
//    double compute_total_generalized_content_with_gradient_and_TGC_projectedHessian(
//            const Eigen::Matrix3Xd &vertices,
//            const Eigen::VectorXi &freeI,
//            const Eigen::Matrix4Xi  &T_free,
//            // output
//            Eigen::VectorXd &generalized_content_list,
//            Eigen::Matrix3Xd &grad,
//            SpMat &Hess) const;

    void set_free_tetI(const Eigen::VectorXi &free_tet_indices) { free_tetI = free_tet_indices; }

private:

    // compute generalized tet volume
    // input:
    //  - v1, v2, v3, v4: four vertices
    double compute_generalized_TetVolume(
            const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
            const Eigen::Vector3d &v3, const Eigen::Vector3d &v4,
            double squared_rest_volume,
            const Eigen::Matrix3d &rest_inverse_EdgeMat
    ) const;

    // compute generalized negative tet volume
    // input:
    //  - v1, v2, v3, v4: four vertices
    double compute_generalized_negative_TetVolume(
            const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
            const Eigen::Vector3d &v3, const Eigen::Vector3d &v4,
            double squared_rest_volume,
            const Eigen::Matrix3d &rest_inverse_EdgeMat
    ) const;

    // compute generalized tet volume with gradient wrt. vert
    // input:
    //  - v1, v2, v3, v4: four vertices
    double compute_generalized_TetVolume_with_gradient(
            const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
            const Eigen::Vector3d &v3, const Eigen::Vector3d &v4,
            double squared_rest_volume,
            const Eigen::Matrix3d &rest_inverse_EdgeMat,
            const Eigen::MatrixXd &pFpx,
            Eigen::Matrix3Xd &grad) const;

    // compute generalized tet volume with gradient and PSD-projected Hessian wrt. vert
    // input:
    //  - vert: four vertices
//    double compute_generalized_TetVolume_with_gradient_projectedHessian(
//            const Eigen::Matrix3Xd &vert,
//            double scaled_squared_rest_volume,
//            double rest_volume,
//            const Eigen::Matrix3d &rest_inverse_EdgeMat,
//            const Eigen::MatrixXd &pFpx,
//            double two_alpha_lambda1,
//            double one_plus_two_alpha_lambda2,
//            double coeff_diag,
//            double coeff_off_diag,
//            Eigen::Matrix3Xd &grad, Eigen::MatrixXd &Hess) const;

    // compute generalized tet volume with gradient and PSD-projected Hessian wrt. vert
    // input:
    //  - v1, v2, v3, v4: four vertices
    double compute_generalized_TetVolume_with_gradient_projected_subtracted_Hessian(
            const Eigen::Vector3d &v1, const Eigen::Vector3d &v2,
            const Eigen::Vector3d &v3, const Eigen::Vector3d &v4,
            double squared_rest_volume, double rest_volume,
            const Eigen::Matrix3d &rest_inverse_EdgeMat,
            const Eigen::MatrixXd &pFpx,
            Eigen::Matrix3Xd &grad, Eigen::MatrixXd &Hess) const;

private:
    // alpha parameter
    double param_alpha;
    // lambda parameters
    double param_lambda1;
    double param_lambda2;
    // k parameter
    double param_k;
    // dimension of target vertices
//    size_t vDim;
    // indices of tetrahedron vertices
    Eigen::Matrix4Xi T;
    // squared edge lengths of auxiliary tetrahedrons
//    Eigen::Matrix6Xd restD;
    // volumes of auxiliary tetrahedrons
    Eigen::VectorXd restA;
    // inverse edge matrices of auxiliary tetrahedrons
    std::vector<Eigen::Matrix3d> rest_invEdgeMat;
    // pFpx: derivative of deformation gradient w.r.t. vertices coordinates
    std::vector<Eigen::MatrixXd> pFpx_list;


    // indices of free tetrahedrons (i.e. tetrahedron with at least one free vertex)
    Eigen::VectorXi free_tetI;

    // squared volumes of auxiliary volumes: (A_rest)^2
    Eigen::VectorXd squared_restA;

    // precomputed coefficients for energy
    // 1 + 2 * alpha * lambda2
    double energy_constant_1;
    // 2 * alpha * lambda1 * k^(2/3)
    double energy_constant_2;
    // alpha * lambda1
    double energy_constant_3;
    // 2 * alpha * lambda2 * k^2 + alpha^2
    double energy_constant_4;

    // precomputed coefficients for gradient
    // 2 * alpha * lambda1
    double grad_constant_1;
    // k^(2/3)
    double k_2_3;

    // coefficients for analytic eigensystem of Hessian
//    double two_alpha_lambda1;
//    double one_plus_two_alpha_lambda2;
//    double coeff_diag;
//    double coeff_off_diag;
//    double coeff_off_diag_subtracted;

};


#endif //TLC_TOTAL_GENERALIZED_CONTENT_3D_H
