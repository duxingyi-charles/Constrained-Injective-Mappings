//
// Created by Charles Du on 3/28/22.
//

#ifndef TLC_STGC_FORMULATION_H
#define TLC_STGC_FORMULATION_H

#include "Total_Generalized_Content.h"

using namespace Eigen;

class sTGC_Formulation {
public:
    // parameter should satisfy: lambda1 + k * lambda2 = 1/2
    sTGC_Formulation(const MatrixXd& rest_vertices, Matrix2Xd init_vertices,
                         Matrix3Xi faces,
                         const VectorXi& handles, const std::string& form,
                         double alpha, double lambda1, double lambda2, double k,
                         bool scale_rest_mesh, bool subtract_total_signed_area,
                         double aspect_ratio_threshold);

    ~sTGC_Formulation() = default;


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
    double compute_energy_with_gradient_projectedHessian(const Eigen::VectorXd &x,
                                                               Eigen::VectorXd &energy_list,
                                                               Eigen::VectorXd &grad, SpMat &Hess);


//    static double computeAlpha(const MatrixXd &restV,
//                               const Matrix2Xd &initV,
//                               const Matrix3Xi &Faces,
//                               const std::string &form, double alphaRatio);

    VectorXd get_x0() const { return x0; }

    std::vector<std::pair<size_t, size_t>> get_boundary_edges() const { return boundary_edges; }

    VectorXi get_freeI() const { return freeI; }

    const Matrix2Xd& get_V() const { return V; }

    const MatrixXd& get_scaled_rest_vertices() const { return scaled_rest_vertices; }

private:
    // x = Flatten(freeV)
    void update_V(const Eigen::VectorXd& x);

private:
    // V indices of triangles
    Matrix3Xi F;
    // indices of free vertices
    VectorXi freeI;
    // indices of free faces (i.e. face with at least one free vertex)
    VectorXi free_faceI;
    // map: V index --> freeV index. If V(i) is not free, indexDict(i) = -1.
    VectorXi indexDict;
    // freeV indices of triangles.
    Matrix3Xi F_free;

    // rest vertices, maybe scaled
    MatrixXd scaled_rest_vertices;

    // boundary edges
    std::vector<std::pair<size_t, size_t>> boundary_edges;

    // TGC energy
    Total_Generalized_Content tgc;

    // current V of target mesh
    Matrix2Xd V;

    // initial variable vector
    VectorXd x0;

    // whether to subtract total signed area
    bool subtract_total_signed_area;

    // whether mesh boundary is fixed
    bool fixed_boundary;
    // signed area of the fixed boundary
    double boundary_signed_area;

    // whether to skip non-free triangles
    bool skip_non_free_triangles;
};


#endif //TLC_STGC_FORMULATION_H
