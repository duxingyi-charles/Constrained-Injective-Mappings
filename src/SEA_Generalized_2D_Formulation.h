//
// Created by Charles Du on 4/27/22.
//

#ifndef TLC_SEA_FORMULATION_H
#define TLC_SEA_FORMULATION_H

#include "Total_Generalized_Content_2D.h"
#include "Arc_Occupancy.h"

#include "timing.h"

using namespace Eigen;

class SEA_Generalized_2D_Formulation {
public:
    // parameter should satisfy: lambda1 + k * lambda2 = 1/2
    SEA_Generalized_2D_Formulation(const MatrixXd& rest_vertices, Matrix2Xd init_vertices,
                     Matrix3Xi faces,
                     const VectorXi& handles, const std::string& form,
                     double alpha, double lambda1, double lambda2, double k,
                     double theta,
                     bool scale_rest_mesh,
                     double aspect_ratio_threshold);

    ~SEA_Generalized_2D_Formulation() = default;


    double compute_energy(const Eigen::VectorXd& x);

    // compute generalized Smooth Excess Area energy,
    // record the decomposition of the energy in energy_list
    // decomposition:
    // - generalized content per triangle
    // - arc segment area per edge
    // - arc occupancy
    double compute_energy(const Eigen::VectorXd& x, Eigen::VectorXd& energy_list);

    double compute_energy_with_gradient(const Eigen::VectorXd& x,
                                        Eigen::VectorXd& grad);


    VectorXd get_x0() const { return x0; }

    std::vector<std::pair<size_t, size_t>> get_boundary_edges() const { return boundary_edges; }

    VectorXi get_freeI() const { return freeI; }

    const Matrix2Xd& get_V() const { return V; }

    const MatrixXd& get_scaled_rest_vertices() const { return scaled_rest_vertices; }

    // whether arc-arc intersection is found in the most recent call of
    // compute_energy() and compute_energy_with_gradient()
    bool has_found_intersection;

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
    Total_Generalized_Content_2D tgc;
    // arc occupancy
    Arc_Occupancy arcOccupancy;

    // current V of target mesh
    Matrix2Xd V;

    // initial variable vector
    VectorXd x0;

    // whether mesh boundary is fixed
    bool fixed_boundary;
    // signed area of the fixed boundary
//    double boundary_signed_area;

    // whether to skip non-free triangles
    bool skip_non_free_triangles;
};


#endif //TLC_SEA_FORMULATION_H
