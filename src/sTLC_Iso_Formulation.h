//
// Created by Charles Du on 2/23/22.
//

#ifndef S_TLC_ISO_FORMULATION_H
#define S_TLC_ISO_FORMULATION_H

#include "Total_Lifted_Content_Isometric.h"

using namespace Eigen;

class sTLC_Iso_Formulation {
public:
    sTLC_Iso_Formulation(const MatrixXd& rest_vertices, Matrix2Xd init_vertices,
        Matrix3Xi faces,
        const VectorXi& handles, const std::string& form,
        double alpha, bool scale_rest_mesh, bool subtract_total_signed_area);

    ~sTLC_Iso_Formulation() = default;


    double compute_energy(const Eigen::VectorXd& x);

    // compute subtracted TLC (isometric) energy,
    // record the decomposition of the energy in energy_list
    // decomposition: subtracted lifted content per triangle
    double compute_energy(const Eigen::VectorXd& x, Eigen::VectorXd& energy_list);

    double compute_energy_with_gradient(const Eigen::VectorXd& x,
        Eigen::VectorXd& grad);

    // PSD-projected Hessian of (TLC_Iso - total signed area)
//    double compute_energy_with_gradient_projectedHessian(const Eigen::VectorXd &x,
//                                                               Eigen::VectorXd &grad, SpMat &Hess);

    // PSD-projected Hessian of (TLC_Iso - total signed area)
//    double compute_energy_with_gradient_projectedHessian(const Eigen::VectorXd &x,
//                                                               Eigen::VectorXd &energy_list,
//                                                               Eigen::VectorXd &grad, SpMat &Hess);


//    static double computeAlpha(const MatrixXd &restV,
//                               const Matrix2Xd &initV,
//                               const Matrix3Xi &Faces,
//                               const std::string &form, double alphaRatio);

    VectorXd get_x0() const { return x0; }

    std::vector<std::pair<size_t, size_t>> get_boundary_edges() const { return boundary_edges; }

    VectorXi get_freeI() const { return freeI; }

    const Matrix2Xd& get_V() const { return V; }

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

    // boundary edges
    std::vector<std::pair<size_t, size_t>> boundary_edges;

    // TLC_Iso energy
    Total_Lifted_Content_Isometric tlc_iso;

    // current V of target mesh
    Matrix2Xd V;

    // initial variable vector
    VectorXd x0;

    // whether to subtract total signed area
    bool subtract_total_signed_area;

};


#endif //S_TLC_ISO_FORMULATION_H
