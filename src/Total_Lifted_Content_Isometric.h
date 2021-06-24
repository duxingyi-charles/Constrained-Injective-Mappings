//
// Created by Charles Du on 6/19/21.
//

#ifndef TLC_TOTAL_LIFTED_CONTENT_ISOMETRIC_H
#define TLC_TOTAL_LIFTED_CONTENT_ISOMETRIC_H

#include "geo_util.h"

#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> SpMat;

class Total_Lifted_Content_Isometric {
public:

    Total_Lifted_Content_Isometric() = default;

    // initialize from rest mesh and parameters: get ready for computing TLC_iso energy and derivatives
    void initialize(const Eigen::MatrixXd &rest_vertices, Eigen::Matrix3Xi faces,
                    const std::string &form, double alpha);

    ~Total_Lifted_Content_Isometric() = default;

    // compute total lifted content (isometric)
    double compute_total_lifted_content_isometric(const Eigen::Matrix2Xd &vertices) const;

    // compute total lifted content (isometric), record lifted content of each triangle
    double compute_total_lifted_content_isometric(const Eigen::Matrix2Xd &vertices,
                                                  Eigen::VectorXd &lifted_content_list) const;

    // compute total lifted content (isometric) and its gradient
    double compute_total_lifted_content_isometric_with_gradient(const Eigen::Matrix2Xd &vertices,
                                                      // output
                                                      Eigen::Matrix2Xd &grad) const;

    // compute total lifted content (isometric) and its gradient,
    // and PSD projected Hessian of (TLC_iso - total signed area) on free vertices
//    double compute_total_lifted_content_with_gradient_and_sTLC_projectedHessian(
//            const Eigen::Matrix2Xd &vertices,
//            const Eigen::VectorXi &freeI,
//            const Eigen::Matrix3Xi  &F_free,
//            // output
//            Eigen::Matrix2Xd &grad,
//            SpMat &Hess) const;

    // compute total lifted content (isometric) and its gradient,
    // and PSD projected Hessian of (TLC_iso - total signed area) on free vertices
//    double compute_total_lifted_content_isometric_with_gradient_and_sTLC_projectedHessian(
//            const Eigen::Matrix2Xd &vertices,
//            const Eigen::VectorXi &freeI,
//            const Eigen::Matrix3Xi  &F_free,
//            // output
//            Eigen::VectorXd &lifted_content_list,
//            Eigen::Matrix2Xd &grad,
//            SpMat &Hess) const;

    void set_free_faceI(const Eigen::VectorXi &free_face_indices) { free_faceI = free_face_indices; }

private:

    // compute lifted triangle area (isometric)
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    double compute_lifted_TriArea_isometric(
            const Eigen::Matrix2Xd &vert,
            const Eigen::Vector3d &Dirichlet_coef,
            double scaled_squared_rest_area) const;

    // compute lifted triangle area (isometric) with gradient wrt. vert
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
    double compute_lifted_TriArea_isometric_with_gradient(
            const Eigen::Matrix2Xd &vert,
            const Eigen::Vector3d &Dirichlet_coef,
            double scaled_squared_rest_area,
            Eigen::Matrix2Xd &grad) const;

    // compute lifted triangle area (isometric) with gradient and Hessian wrt. vert
    // input:
    //  - vert: three vertices
    //  - r: squared edge lengths of aux triangle
//    static double compute_lifted_TriArea_isometric_with_gradient_Hessian(
//            const Eigen::Matrix2Xd &vert, const Eigen::Vector3d &r,
//            Eigen::Matrix2Xd &grad, Eigen::MatrixXd &Hess);

private:
    // alpha parameter
    double param_alpha;
    // dimension of target vertices
//    size_t vDim;
    // indices of triangle vertices
    Eigen::Matrix3Xi F;
    // squared edge lengths of auxiliary triangles
    Eigen::Matrix3Xd restD;
    // areas of auxiliary triangles
    Eigen::VectorXd restA;

    // indices of free faces (i.e. face with at least one free vertex)
    Eigen::VectorXi free_faceI;

    // scaled squared areas of auxiliary triangles: (alpha/2 + alpha^2)*(A_rest)^2
    Eigen::VectorXd scaled_squared_restA;
    // scaled Dirichlet coefficients
    Eigen::Matrix3Xd scaled_Dirichlet_coef;

};


#endif //TLC_TOTAL_LIFTED_CONTENT_ISOMETRIC_H
