//
// Created by Charles Du on 3/29/22.
//

#ifndef TLC_DEFORMATION_GRADIENT_UTIL_H
#define TLC_DEFORMATION_GRADIENT_UTIL_H

#include <Eigen/Core>


// compute edge matrix for a triangle
void compute_edge_matrix(const Eigen::VectorXd& v1,
                         const Eigen::VectorXd& v2,
                         const Eigen::VectorXd& v3,
                         Eigen::Matrix2d& edgeMat);


// compute edge matrix for a tetrahedron
void compute_edge_matrix(const Eigen::Vector3d& v1,
                         const Eigen::Vector3d& v2,
                         const Eigen::Vector3d& v3,
                         const Eigen::Vector3d& v4,
                         Eigen::Matrix3d& edgeMat);

// deformation gradient
// = edge_matrix(target simplex) * Inverse(edge_matrix(auxiliary simplex))

// vectorize a matrix
void vectorize(const Eigen::MatrixXd& A,
               Eigen::VectorXd& vecA);

// compute flattened 4th-order tensor pFpx for a triangle
// input: inverse of auxiliary triangle's edge matrix
// output: derivative of deformation gradient w.r.t. target vertices,
// flattened to a matrix
void compute_flatten_pFpx(const Eigen::Matrix2d& invEdgeMat,
                          Eigen::MatrixXd& pFpx);

// compute flattened 4th-order tensor pFpx for a tetrahedron
// input: inverse of auxiliary tetrahedron's edge matrix
// output: derivative of deformation gradient w.r.t. target vertices,
// flattened to a matrix
void compute_flatten_pFpx(const Eigen::Matrix3d& invEdgeMat,
                          Eigen::MatrixXd& pFpx);

// compute rotation-variant SVD for 2x2 matrix
void compute_SVD(const Eigen::Matrix2d& mat,
                 Eigen::Matrix2d& U,
                 Eigen::Vector2d& singular_values,
                 Eigen::Matrix2d& V);

// compute rotation-variant SVD for 3x3 matrix
void compute_SVD(const Eigen::Matrix3d& mat,
                 Eigen::Matrix3d& U,
                 Eigen::Vector3d& singular_values,
                 Eigen::Matrix3d& V);

#endif //TLC_DEFORMATION_GRADIENT_UTIL_H
