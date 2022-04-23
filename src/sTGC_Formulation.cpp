//
// Created by Charles Du on 3/28/22.
//

#include "sTGC_Formulation.h"
#include <iostream>

sTGC_Formulation::sTGC_Formulation(const MatrixXd &rest_vertices, Matrix2Xd init_vertices, Matrix3Xi faces,
                                   const VectorXi &handles, const std::string &form, double alpha, double lambda1,
                                   double lambda2, double k, bool scale_rest_mesh, bool subtract_total_signed_area) :
        F(std::move(faces)), V(std::move(init_vertices)), subtract_total_signed_area(subtract_total_signed_area),
        skip_non_free_triangles(false)
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
    double init_area = compute_total_signed_mesh_area(V,F);
    double rest_area = compute_total_unsigned_area(rest_vertices,F);
    double rest_scale = 1;
    if (scale_rest_mesh && handles.size() > 1 && init_area > 0) {
        rest_scale = sqrt(init_area / rest_area);
        std::cout << "scale rest mesh by " << rest_scale << std::endl;
    }

    // print formulation parameters
    std::cout << "form: " << form << std::endl;
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "lambda1: " << lambda1 << std::endl;
    std::cout << "lambda2: " << lambda2 << std::endl;
    std::cout << "k: " << k << std::endl;

    // compute x0 from initV
    x0.resize(vDim * freeI.size());
    for (auto i = 0; i < freeI.size(); ++i) {
        int vi = freeI(i);
        for (int j = 0; j < vDim; ++j) {
            x0(i * vDim + j) = V(j, vi);
        }
    }

    // initialize TGC
    scaled_rest_vertices = rest_scale * rest_vertices;
    tgc.initialize(scaled_rest_vertices, F, form, alpha, lambda1, lambda2, k);
//    std::cout << "free_faceI: " << std::endl;
//    std::cout << free_faceI << std::endl;
    if (skip_non_free_triangles) {
        tgc.set_free_faceI(free_faceI);
    }

    // extract boundary edges
    extract_mesh_boundary_edges(F, boundary_edges);

    // check if the boundary is fixed
    fixed_boundary = true;
    for (const auto& edge: boundary_edges) {
        if (freeQ[edge.first] || freeQ[edge.second]) {
            fixed_boundary = false;
            break;
        }
    }
    boundary_signed_area = 0;
    if (fixed_boundary) {
        boundary_signed_area = init_area;
    }
}

double sTGC_Formulation::compute_energy(const VectorXd &x) {
    update_V(x);

    // TGC
    double tgc_energy = tgc.compute_total_generalized_content(V);

    if (subtract_total_signed_area) {
        if (fixed_boundary) {
            return tgc_energy - boundary_signed_area;
        }
        // convert V to a vector of Point
        std::vector<Point> vertices(V.cols());
        for (int i = 0; i < V.cols(); ++i) {
            vertices[i] = V.col(i);
        }
        double total_signed_area = compute_total_signed_area(vertices, boundary_edges);
        //
        return tgc_energy - total_signed_area;
    }
    else {
        return tgc_energy;
    }
}


void sTGC_Formulation::update_V(const VectorXd &x) {
    int vDim = 2;
    for (auto i = 0; i < freeI.size(); ++i) {
        for (auto j = 0; j < vDim; ++j) {
            V(j, freeI(i)) = x[i * vDim + j];
        }
    }
}

double sTGC_Formulation::compute_energy(const VectorXd &x, VectorXd &energy_list) {
    update_V(x);

    if (subtract_total_signed_area) { // sTGC
        return tgc.compute_total_generalized_negative_content(V, energy_list);
    }
    else { // TGC
        return tgc.compute_total_generalized_content(V, energy_list);
    }
}

double sTGC_Formulation::compute_energy_with_gradient(const VectorXd &x, VectorXd &grad) {
    update_V(x);

    // TGC
    Matrix2Xd tgc_grad;
    double tgc_energy = tgc.compute_total_generalized_content_with_gradient(V, tgc_grad);

    if (subtract_total_signed_area && !fixed_boundary) {
        //
        std::vector<Point> vertices(V.cols());
        for (int i = 0; i < V.cols(); ++i) {
            vertices[i] = V.col(i);
        }
        Matrix2Xd total_signed_area_grad;
        total_signed_area_grad.setZero(2, vertices.size());
        double total_signed_area = compute_total_signed_area_with_gradient(
                vertices, boundary_edges, total_signed_area_grad);

        // full gradient
        Matrix2Xd m_grad = tgc_grad - total_signed_area_grad;

        // gradient of free vertices
        grad.resize(freeI.size() * m_grad.rows());
        for (int i = 0; i < freeI.size(); ++i) {
            grad(2*i) = m_grad(0,freeI(i));
            grad(2*i+1) = m_grad(1,freeI(i));
        }

        // energy
        return tgc_energy - total_signed_area;
    }
    else {
        // gradient of free vertices
        grad.resize(freeI.size() * tgc_grad.rows());
        for (int i = 0; i < freeI.size(); ++i) {
            grad(2*i) = tgc_grad(0,freeI(i));
            grad(2*i+1) = tgc_grad(1,freeI(i));
        }

        // energy
        if (subtract_total_signed_area) { // the boundary must be fixed
            return tgc_energy - boundary_signed_area;
        } else {
            return tgc_energy;
        }
    }
}

double sTGC_Formulation::compute_energy_with_gradient_projectedHessian(const VectorXd &x, VectorXd &energy_list,
                                                                       VectorXd &grad, SpMat &Hess) {
    update_V(x);

    // TGC energy, TGC full gradient,
    // and PSD projected Hessian of (TGC - total signed area) on free vertices
    VectorXd generalized_content_list;
    Matrix2Xd tgc_grad;
    double tgc_energy;
    if (!subtract_total_signed_area) {
        tgc_energy = tgc.compute_total_generalized_content_with_gradient_and_TGC_projectedHessian(V, freeI,
                                                                                                          F_free,
                                                                                                          generalized_content_list,
                                                                                                          tgc_grad,
                                                                                                          Hess);
        energy_list = generalized_content_list;

    } else {
        tgc_energy = tgc.compute_total_generalized_content_with_gradient_and_sTGC_projectedHessian(V, freeI,
                                                                                                   F_free,
                                                                                                   generalized_content_list,
                                                                                                   tgc_grad,
                                                                                                   Hess);
        // signed area
        VectorXd signed_area_list;
        compute_signed_tri_areas(V, F, signed_area_list);

        // fill the energy decomposition into energy_list
        energy_list.resize(generalized_content_list.size());
        for (int i = 0; i < generalized_content_list.size(); ++i) {
            energy_list(i) = generalized_content_list(i) - signed_area_list(i);
        }
    }

    if (subtract_total_signed_area && !fixed_boundary) {
        //
        std::vector<Point> vertices(V.cols());
        for (int i = 0; i < V.cols(); ++i) {
            vertices[i] = V.col(i);
        }
        Matrix2Xd total_signed_area_grad;
        total_signed_area_grad.setZero(2, vertices.size());
        double total_signed_area = compute_total_signed_area_with_gradient(
                vertices, boundary_edges, total_signed_area_grad);

        // full gradient
        Matrix2Xd m_grad = tgc_grad - total_signed_area_grad;

        // gradient of free vertices
        grad.resize(freeI.size() * m_grad.rows());
        for (int i = 0; i < freeI.size(); ++i) {
            grad(2*i) = m_grad(0,freeI(i));
            grad(2*i+1) = m_grad(1,freeI(i));
        }

        // energy
        return tgc_energy - total_signed_area;
    }
    else {
        // gradient of free vertices
        grad.resize(freeI.size() * tgc_grad.rows());
        for (int i = 0; i < freeI.size(); ++i) {
            grad(2*i) = tgc_grad(0,freeI(i));
            grad(2*i+1) = tgc_grad(1,freeI(i));
        }

        // energy
        if (subtract_total_signed_area) { // the boundary must be fixed
            return tgc_energy - boundary_signed_area;
        } else {
            return tgc_energy;
        }
    }
}
