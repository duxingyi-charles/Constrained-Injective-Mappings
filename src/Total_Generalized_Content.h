//
// Created by Charles Du on 3/28/22.
//

#ifndef TLC_TOTAL_GENERALIZED_CONTENT_H
#define TLC_TOTAL_GENERALIZED_CONTENT_H

#include "geo_util.h"

#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> SpMat;

// Total Smooth Generalized Content in 2D
class Total_Generalized_Content {
public:
    Total_Generalized_Content() = default;

    // initialize from rest mesh and parameters: get ready for computing TGC energy and derivatives
    // need to satisfy lambda1 + k * lambda2 = 1/2
    void initialize(const Eigen::MatrixXd &rest_vertices, const Eigen::Matrix3Xi &faces,
                    const std::string &form, double alpha, double lambda1, double lambda2, double k);

    ~Total_Generalized_Content() = default;

    // compute total generalized content
    double compute_total_generalized_content(const Eigen::Matrix2Xd &vertices) const;

    // compute total generalized content, record generalized content of each triangle
    double compute_total_generalized_content(const Eigen::Matrix2Xd &vertices,
                                                  Eigen::VectorXd &generalized_content_list) const;

    // compute generalized negative content, record (generalized content - signed content) of each triangle
    double compute_total_generalized_negative_content(const Eigen::Matrix2Xd &vertices,
                                                      Eigen::VectorXd &generalized_neg_content_list) const;

    // compute total generalized content and its gradient
    double compute_total_generalized_content_with_gradient(const Eigen::Matrix2Xd &vertices,
                                                                Eigen::Matrix2Xd &grad) const;

    // compute total generalized content and its gradient,
    // and PSD projected Hessian of (TGC - total signed area) on free vertices
//    double compute_total_generalized_content_with_gradient_and_sTGC_projectedHessian(
//            const Eigen::Matrix2Xd &vertices,
//            const Eigen::VectorXi &freeI,
//            const Eigen::Matrix3Xi  &F_free,
//            // output
//            Eigen::Matrix2Xd &grad,
//            SpMat &Hess) const;

    // compute total generalized content and its gradient,
    // and PSD projected Hessian of (TGC - total signed area) on free vertices
    double compute_total_generalized_content_with_gradient_and_sTGC_projectedHessian(
            const Eigen::Matrix2Xd &vertices,
            const Eigen::VectorXi &freeI,
            const Eigen::Matrix3Xi  &F_free,
            // output
            Eigen::VectorXd &generalized_content_list,
            Eigen::Matrix2Xd &grad,
            SpMat &Hess) const;

    // compute total generalized content and its gradient,
    // and PSD projected Hessian of TGC on free vertices
    double compute_total_generalized_content_with_gradient_and_TGC_projectedHessian(
            const Eigen::Matrix2Xd &vertices,
            const Eigen::VectorXi &freeI,
            const Eigen::Matrix3Xi  &F_free,
            // output
            Eigen::VectorXd &generalized_content_list,
            Eigen::Matrix2Xd &grad,
            SpMat &Hess) const;

    void set_free_faceI(const Eigen::VectorXi &free_face_indices) { free_faceI = free_face_indices; }

private:

    // compute generalized triangle area
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    double compute_generalized_TriArea(
            const Eigen::Vector2d &v1,
            const Eigen::Vector2d &v2,
            const Eigen::Vector2d &v3,
            const Eigen::Vector3d &Dirichlet_coef,
            double scaled_squared_rest_area) const;

    // compute generalized negative triangle area
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    double compute_generalized_negative_TriArea(
            const Eigen::Vector2d &v1,
            const Eigen::Vector2d &v2,
            const Eigen::Vector2d &v3,
            const Eigen::Vector3d &Dirichlet_coef,
            double scaled_squared_rest_area) const;

    // compute generalized triangle area with gradient wrt. vert
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    double compute_generalized_TriArea_with_gradient(
            const Eigen::Matrix2Xd &vert,
            const Eigen::Vector3d &Dirichlet_coef,
            double scaled_squared_rest_area,
            Eigen::Matrix2Xd &grad) const;

    // compute generalized triangle area with gradient and PSD-projected Hessian wrt. vert
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    double compute_generalized_TriArea_with_gradient_projectedHessian(
            const Eigen::Matrix2Xd &vert,
            const Eigen::Vector3d &Dirichlet_coef,
            double scaled_squared_rest_area,
            double rest_area,
            const Eigen::Matrix2d &rest_inverse_EdgeMat,
            const Eigen::MatrixXd &pFpx,
            double two_alpha_lambda1,
            double one_plus_two_alpha_lambda2,
            double coeff_diag,
            double coeff_off_diag,
            Eigen::Matrix2Xd &grad, Eigen::MatrixXd &Hess) const;

    // compute generalized triangle area with gradient and PSD-projected Hessian wrt. vert
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    double compute_generalized_TriArea_with_gradient_projected_subtracted_Hessian(
            const Eigen::Vector2d &v1, const Eigen::Vector2d &v2, const Eigen::Vector2d &v3,
            const Eigen::Vector3d &Dirichlet_coef,
            double scaled_squared_rest_area,
            double rest_area,
            const Eigen::Matrix2d &rest_inverse_EdgeMat,
            const Eigen::MatrixXd &pFpx,
            double two_alpha_lambda1,
            double one_plus_two_alpha_lambda2,
            double coeff_diag,
            double coeff_off_diag,
            Eigen::Matrix2Xd &grad, Eigen::MatrixXd &Hess) const;

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
    // indices of triangle vertices
    Eigen::Matrix3Xi F;
    // squared edge lengths of auxiliary triangles
    Eigen::Matrix3Xd restD;
    // areas of auxiliary triangles
    Eigen::VectorXd restA;
    // inverse edge matrices of auxiliary triangles
    std::vector<Eigen::Matrix2d> rest_invEdgeMat;
    // pFpx: derivative of deformation gradient w.r.t. vertices coordinates
    std::vector<Eigen::MatrixXd> pFpx_list;


    // indices of free faces (i.e. face with at least one free vertex)
    Eigen::VectorXi free_faceI;

    // coefficient for the squared areas of target triangles: (1+2*alpha*lambda2)
    double squared_targetA_scale_coeff;
    // scaled squared areas of auxiliary triangles: (2*alpha*lambda2*k^2 + alpha^2)*(A_rest)^2
    Eigen::VectorXd scaled_squared_restA;
    // scaled Dirichlet coefficients: (alpha*lambda1/4)(d_{i+1} + d{i+2} - d_i)
    // where d_i is the squared length of i-th edge of the auxiliary triangle
    Eigen::Matrix3Xd scaled_Dirichlet_coef;

    // coefficients for analytic eigensystem of Hessian
    double two_alpha_lambda1;
    double one_plus_two_alpha_lambda2;
    double coeff_diag;
    double coeff_off_diag;
    double coeff_off_diag_subtracted;
};


#endif //TLC_TOTAL_GENERALIZED_CONTENT_H
