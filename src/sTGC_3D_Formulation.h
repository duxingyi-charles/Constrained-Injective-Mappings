//
// Created by Charles Du on 4/5/22.
//

#ifndef TLC_STGC_3D_FORMULATION_H
#define TLC_STGC_3D_FORMULATION_H

#include "Total_Generalized_Content_3D.h"
using namespace Eigen;


class sTGC_3D_Formulation {
public:
    // parameter should satisfy: 6 * k^(1/3) * lambda1 + 2 * k * lambda2 = 1
    sTGC_3D_Formulation(const Matrix3Xd& rest_vertices, Matrix3Xd init_vertices,
                     Matrix4Xi tets,
                     const VectorXi& handles, const std::string& form,
                     double alpha, double lambda1, double lambda2, double k,
                     bool scale_rest_mesh, bool subtract_total_signed_content);

    ~sTGC_3D_Formulation() = default;


    double compute_energy(const Eigen::VectorXd& x);

    // compute subtracted TGC energy,
    // record the decomposition of the energy in energy_list
    // decomposition: subtracted generalized content per triangle
    double compute_energy(const Eigen::VectorXd& x, Eigen::VectorXd& energy_list);

    double compute_energy_with_gradient(const Eigen::VectorXd& x,
                                        Eigen::VectorXd& grad);

    // PSD-projected Hessian of (TGC - total signed area)
//    double compute_energy_with_gradient_projectedHessian(const Eigen::VectorXd &x,
//                                                               Eigen::VectorXd &grad, SpMat &Hess);

    // PSD-projected Hessian of (TGC - total signed area)
//    double compute_energy_with_gradient_projectedHessian(const Eigen::VectorXd &x,
//                                                         Eigen::VectorXd &energy_list,
//                                                         Eigen::VectorXd &grad, SpMat &Hess);


//    static double computeAlpha(const MatrixXd &restV,
//                               const Matrix2Xd &initV,
//                               const Matrix3Xi &Faces,
//                               const std::string &form, double alphaRatio);

    VectorXd get_x0() const { return x0; }

    std::vector<std::array<int,3>> get_boundary_triangles() const { return boundary_triangles; }

    VectorXi get_freeI() const { return freeI; }

    const Matrix3Xd& get_V() const { return V; }

private:
    // x = Flatten(freeV)
    void update_V(const Eigen::VectorXd& x);

private:
    // V indices of tetrahedrons
    Matrix4Xi T;
    // indices of free vertices
    VectorXi freeI;
    // indices of free tetrahedrons (i.e. tetrahedron with at least one free vertex)
    VectorXi free_tetI;
    // map: V index --> freeV index. If V(i) is not free, indexDict(i) = -1.
    VectorXi indexDict;
    // freeV indices of tetrahedrons.
    Matrix4Xi T_free;

    // boundary triangles
    std::vector<std::array<int,3>> boundary_triangles;

    // TGC energy
    Total_Generalized_Content_3D tgc;

    // current V of target mesh
    Matrix3Xd V;

    // initial variable vector
    VectorXd x0;

    // whether to subtract total signed volume
    bool subtract_total_signed_content;

    // whether mesh boundary is fixed
    bool fixed_boundary;
    // signed content of the fixed boundary
    double boundary_signed_content;

    // whether to skip non-free tetrahedrons
    bool skip_non_free_tets;
};


#endif //TLC_STGC_3D_FORMULATION_H
