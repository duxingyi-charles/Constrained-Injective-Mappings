//
// Created by Charles Du on 6/19/21.
//

#include "Arc_Overlap_Iso_Formulation.h"
#include <iostream>

Arc_Overlap_Iso_Formulation::Arc_Overlap_Iso_Formulation(const MatrixXd &rest_vertices, Matrix2Xd init_vertices,
                                                         Matrix3Xi faces, const VectorXi &handles,
                                                         const std::string &form, double alpha, double theta,
                                                         bool scale_rest_mesh) :
        F(std::move(faces)), V(std::move(init_vertices)),
        arcOccupancy(theta),
        // we assume there are intersection at the beginning
        has_found_intersection(true)
{
    // compute freeI: indices of free vertices
    int nV = rest_vertices.cols();
    int vDim = 2;

    std::vector<bool> freeQ(nV, true);
    for (auto i = 0; i < handles.size(); ++i) {
        freeQ[handles(i)] = false;
    }
    freeI.resize(nV - handles.size());
    int ii = 0;
    for (int i = 0; i < nV; ++i) {
        if (freeQ[i]) {
            freeI[ii] = i;
            ++ii;
        }
    }
    std::sort(freeI.data(), freeI.data() + freeI.size());

    // experimental: compute free face Indices
    // free face: a triangle with at least one free vertices
    int nF = F.cols();
    int n_free_face = 0;
    std::vector<bool> free_faceQ(nF, false);
    for (int i = 0; i < nF; ++i) {
        int i1 = F(0,i);
        int i2 = F(1,i);
        int i3 = F(2,i);
        if (freeQ[i1] || freeQ[i2] || freeQ[i3]) {
            free_faceQ[i] = true;
            ++n_free_face;
        }
    }
//    std::cout << "n_free_face = " << n_free_face << std::endl;
    free_faceI.resize(n_free_face);
    ii = 0;
    for (int i = 0; i < nF; ++i) {
        if (free_faceQ[i]) {
            free_faceI[ii] = i;
            ++ii;
        }
    }


    // compute indexDict and F_free
    indexDict = VectorXi::Constant(nV, -1);
    for (auto i = 0; i < freeI.size(); ++i) {
        indexDict(freeI(i)) = i;
    }

    F_free.resize(F.rows(), F.cols());
    for (auto i = 0; i < F.cols(); ++i) {
        for (auto j = 0; j < F.rows(); ++j) {
            F_free(j, i) = indexDict(F(j, i));
        }
    }

    // (optional) scale rest mesh to have the same area as init mesh
    double init_area = compute_total_signed_mesh_area(init_vertices,faces);
    double rest_area = compute_total_unsigned_area(rest_vertices,faces);
    double rest_scale = 1;
    if (scale_rest_mesh && handles.size() > 1 && init_area > 0) {
        rest_scale = sqrt(init_area / rest_area);
        std::cout << "scale rest mesh by " << rest_scale << std::endl;
    }

    // compute alpha if it's not explicitly given
    std::cout << "form: " << form << std::endl;
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "theta: " << theta << std::endl;

    // compute x0 from initV
    x0.resize(vDim * freeI.size());
    for (auto i = 0; i < freeI.size(); ++i) {
        int vi = freeI(i);
        for (int j = 0; j < vDim; ++j) {
            x0(i * vDim + j) = V(j, vi);
        }
    }

    // initialize TLC_Iso
    tlc_iso.initialize(rest_scale * rest_vertices, F, form, alpha);
//    std::cout << "free_faceI: " << std::endl;
//    std::cout << free_faceI << std::endl;
    tlc_iso.set_free_faceI(free_faceI);

    // extract boundary edges
    extract_mesh_boundary_edges(F, boundary_edges);
}

double Arc_Overlap_Iso_Formulation::compute_energy(const VectorXd &x) {
    update_V(x);

    // TLC (isometric)
    double tlc_iso_energy = tlc_iso.compute_total_lifted_content_isometric(V);

    //
    std::vector<Point> vertices(V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        vertices[i] = V.col(i);
    }
    double arc_segment_area = arcOccupancy.compute_arc_curve_segment_area(vertices, boundary_edges);
    double arc_occupancy = arcOccupancy.compute_arc_occupancy(vertices, boundary_edges);
    has_found_intersection = arcOccupancy.has_found_intersection;

    return tlc_iso_energy + arc_segment_area - arc_occupancy;
}

void Arc_Overlap_Iso_Formulation::update_V(const VectorXd &x) {
    int vDim = 2;
    for (auto i = 0; i < freeI.size(); ++i) {
        for (auto j = 0; j < vDim; ++j) {
            V(j, freeI(i)) = x[i * vDim + j];
        }
    }
}

double Arc_Overlap_Iso_Formulation::compute_energy(const VectorXd &x, VectorXd &energy_list) {
    update_V(x);

    // TLC
    VectorXd lifted_content_list;
    double tlc_iso_energy = tlc_iso.compute_total_lifted_content_isometric(V,lifted_content_list);

    //
    std::vector<Point> vertices(V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        vertices[i] = V.col(i);
    }
    VectorXd arc_curve_segment_area_list;
    double arc_segment_area = arcOccupancy.compute_arc_curve_segment_area(vertices, boundary_edges,
                                                                          arc_curve_segment_area_list);
    double arc_occupancy = arcOccupancy.compute_arc_occupancy(vertices, boundary_edges);
    has_found_intersection = arcOccupancy.has_found_intersection;

    // fill the energy decomposition into energy_list
    energy_list.resize(lifted_content_list.size()+arc_curve_segment_area_list.size()+1);
    int ii = 0;
    for (int i = 0; i < lifted_content_list.size(); ++i) {
        energy_list(ii) = lifted_content_list(i);
        ++ii;
    }
    for (int i = 0; i < arc_curve_segment_area_list.size(); ++i) {
        energy_list(ii) = arc_curve_segment_area_list(i);
        ++ii;
    }
    energy_list(ii) = -arc_occupancy;


    return tlc_iso_energy + arc_segment_area - arc_occupancy;
}

double Arc_Overlap_Iso_Formulation::compute_energy_with_gradient(const VectorXd &x, VectorXd &grad) {
    update_V(x);

    // TLC (isometric)
    Matrix2Xd tlc_grad;
    double tlc_energy = tlc_iso.compute_total_lifted_content_isometric_with_gradient(V, tlc_grad);

    //
    std::vector<Point> vertices(V.cols());
    for (int i = 0; i < V.cols(); ++i) {
        vertices[i] = V.col(i);
    }
    Matrix2Xd arc_seg_grad;
    double arc_segment_area =arcOccupancy.compute_arc_curve_segment_area_with_gradient(vertices, boundary_edges,
                                                                                       arc_seg_grad);
    Matrix2Xd arc_occupancy_grad;
    double arc_occupancy = arcOccupancy.compute_arc_occupancy_with_gradient(vertices, boundary_edges,
                                                                            arc_occupancy_grad);
    has_found_intersection = arcOccupancy.has_found_intersection;

    // full gradient
    Matrix2Xd m_grad = tlc_grad + arc_seg_grad - arc_occupancy_grad;

    // gradient of free vertices
    grad.resize(freeI.size() * m_grad.rows());
    for (int i = 0; i < freeI.size(); ++i) {
        grad(2*i) = m_grad(0,freeI(i));
        grad(2*i+1) = m_grad(1,freeI(i));
    }

    // energy
    return tlc_energy + arc_segment_area - arc_occupancy;
}
