//
// Created by Charles Du on 3/29/22.
//

#include "deformation_gradient_util.h"
#include <Eigen/SVD>
#include <Eigen/LU> // function determinant()

void compute_edge_matrix(const Eigen::VectorXd &v1, const Eigen::VectorXd &v2, const Eigen::VectorXd &v3,
                         Eigen::Matrix2d &edgeMat) {
    if (v1.rows() == 2) { // embedded in 2D
        edgeMat.col(0) = v2 - v1;
        edgeMat.col(1) = v3 - v1;
    }
    else { // embedded in 3D
        double x1 = (v2-v1).norm();
        double x2 = (v2-v1).dot(v3-v1)/x1;
        double y2 = sqrt((v3-v1).squaredNorm() - x2*x2);
        edgeMat << x1, x2,
                    0, y2;
    }
}

void compute_edge_matrix(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d &v3,
                         const Eigen::Vector3d &v4, Eigen::Matrix3d &edgeMat) {
    edgeMat.col(0) = v2 - v1;
    edgeMat.col(1) = v3 - v1;
    edgeMat.col(2) = v4 - v1;
}

void vectorize(const Eigen::MatrixXd &A, Eigen::VectorXd &vecA) {
    vecA.resize(A.cols() * A.rows());
    int index = 0;
    for (int y=0; y<A.cols(); ++y) {
        for (int x = 0; x < A.rows(); ++x, ++index) {
            vecA[index] = A(x,y);
        }
    }
}

void vectorize2x2(const Eigen::Matrix2d &A, Eigen::Vector4d &vecA) {
    vecA[0] = A(0,0);
    vecA[1] = A(1,0);
    vecA[2] = A(0,1);
    vecA[3] = A(1,1);
}

void vectorize3x3(const Eigen::Matrix3d& A, Vector9d& vecA) {
    vecA[0] = A(0,0);
    vecA[1] = A(1,0);
    vecA[2] = A(2,0);
    vecA[3] = A(0,1);
    vecA[4] = A(1,1);
    vecA[5] = A(2,1);
    vecA[6] = A(0,2);
    vecA[7] = A(1,2);
    vecA[8] = A(2,2);
}


void compute_flatten_pFpx(const Eigen::Matrix2d &invEdgeMat,
                          Eigen::MatrixXd& pFpx) {
    double a = invEdgeMat(0,0);
    double b = invEdgeMat(0,1);
    double c = invEdgeMat(1,0);
    double d = invEdgeMat(1,1);
    double ac = -a - c;
    double bd = -b - d;
    // F is 2*2 matrix, x is 3*2 matrix
    pFpx.setZero(4, 6);
    pFpx(0,0) = ac;
    pFpx(2,0) = bd;
    pFpx(1,1) = ac;
    pFpx(3,1) = bd;
    pFpx(0,2) = a;
    pFpx(2,2) = b;
    pFpx(1,3) = a;
    pFpx(3,3) = b;
    pFpx(0,4) = c;
    pFpx(2,4) = d;
    pFpx(1,5) = c;
    pFpx(3,5) = d;
}

void compute_flatten_pFpx(const Eigen::Matrix3d &invEdgeMat, Eigen::MatrixXd &pFpx) {
    double a = invEdgeMat(0, 0);
    double b = invEdgeMat(0, 1);
    double c = invEdgeMat(0, 2);
    double d = invEdgeMat(1, 0);
    double e = invEdgeMat(1, 1);
    double f = invEdgeMat(1, 2);
    double g = invEdgeMat(2, 0);
    double h = invEdgeMat(2, 1);
    double i = invEdgeMat(2, 2);
    double adg = -a - d - g;
    double beh = -b - e - h;
    double cfi = -c - f - i;
    // F is 3*3 matrix, x is 4*3 matrix
    pFpx.setZero(9, 12);
    pFpx(0, 0) = adg;
    pFpx(3, 0) = beh;
    pFpx(6, 0) = cfi;
    pFpx(1, 1) = adg;
    pFpx(4, 1) = beh;
    pFpx(7, 1) = cfi;
    pFpx(2, 2) = adg;
    pFpx(5, 2) = beh;
    pFpx(8, 2) = cfi;

    pFpx(0, 3) = a;
    pFpx(3, 3) = b;
    pFpx(6, 3) = c;
    pFpx(1, 4) = a;
    pFpx(4, 4) = b;
    pFpx(7, 4) = c;
    pFpx(2, 5) = a;
    pFpx(5, 5) = b;
    pFpx(8, 5) = c;

    pFpx(0, 6) = d;
    pFpx(3, 6) = e;
    pFpx(6, 6) = f;
    pFpx(1, 7) = d;
    pFpx(4, 7) = e;
    pFpx(7, 7) = f;
    pFpx(2, 8) = d;
    pFpx(5, 8) = e;
    pFpx(8, 8) = f;

    pFpx(0, 9) = g;
    pFpx(3, 9) = h;
    pFpx(6, 9) = i;
    pFpx(1, 10) = g;
    pFpx(4, 10) = h;
    pFpx(7, 10) = i;
    pFpx(2, 11) = g;
    pFpx(5, 11) = h;
    pFpx(8, 11) = i;
}

void compute_SVD(const Eigen::Matrix2d &mat, Eigen::Matrix2d &U, Eigen::Vector2d &singular_values, Eigen::Matrix2d &V) {
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    singular_values = svd.singularValues();
    U = svd.matrixU();
    V = svd.matrixV();
    double detU = U.determinant();
    double detV = V.determinant();
    if (detU < 0 && detV > 0) {
        U.col(1) *= -1;
        singular_values[1] *= -1;
    }
    else if (detU > 0 && detV < 0) {
        V.col(1) *= -1;
        singular_values[1] *= -1;
    }
}

void compute_SVD(const Eigen::Matrix3d &mat, Eigen::Matrix3d &U, Eigen::Vector3d &singular_values, Eigen::Matrix3d &V) {
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    singular_values = svd.singularValues();
    U = svd.matrixU();
    V = svd.matrixV();
    double detU = U.determinant();
    double detV = V.determinant();
    if (detU < 0 && detV > 0) {
        U.col(2) *= -1;
        singular_values[2] *= -1;
    }
    else if (detU > 0 && detV < 0) {
        V.col(2) *= -1;
        singular_values[2] *= -1;
    }
}

void compute_tri_mesh_singular_values(const Eigen::MatrixXd &rest_vertices, const Eigen::Matrix2Xd &target_vertices,
                                      const Eigen::Matrix3Xi &F, Eigen::Matrix2Xd &singular_values) {
    Eigen::Matrix2d rest_edgeMat;
    Eigen::Matrix2d target_edgeMat;
    Eigen::Matrix2d deformation_gradient;
    Eigen::Matrix2d U, V;
    Eigen::Vector2d s;

    singular_values.resize(2, F.cols());
    for (int i = 0; i < F.cols(); ++i) {
        int i1 = F(0,i);
        int i2 = F(1,i);
        int i3 = F(2,i);

        compute_edge_matrix(rest_vertices.col(i1), rest_vertices.col(i2), rest_vertices.col(i3),
                            rest_edgeMat);
        compute_edge_matrix(target_vertices.col(i1), target_vertices.col(i2), target_vertices.col(i3),
                            target_edgeMat);
        deformation_gradient = target_edgeMat * rest_edgeMat.inverse();
        compute_SVD(deformation_gradient, U, s, V);
        singular_values.col(i) = s;
    }
}

void compute_tet_mesh_singular_values(const Eigen::Matrix3Xd &rest_vertices, const Eigen::Matrix3Xd &target_vertices,
                                      const Eigen::Matrix4Xi &Tets, Eigen::Matrix3Xd &singular_values) {
    Eigen::Matrix3d rest_edgeMat;
    Eigen::Matrix3d target_edgeMat;
    Eigen::Matrix3d deformation_gradient;
    Eigen::Matrix3d U, V;
    Eigen::Vector3d s;

    singular_values.resize(3, Tets.cols());
    for (int i = 0; i < Tets.cols(); ++i) {
        int i1 = Tets(0,i);
        int i2 = Tets(1,i);
        int i3 = Tets(2,i);
        int i4 = Tets(3,i);

        compute_edge_matrix(rest_vertices.col(i1), rest_vertices.col(i2), rest_vertices.col(i3),
                            rest_vertices.col(i4),rest_edgeMat);
        compute_edge_matrix(target_vertices.col(i1), target_vertices.col(i2), target_vertices.col(i3),
                            target_vertices.col(i4),target_edgeMat);
        deformation_gradient = target_edgeMat * rest_edgeMat.inverse();
        compute_SVD(deformation_gradient, U, s, V);
        singular_values.col(i) = s;
    }

}



