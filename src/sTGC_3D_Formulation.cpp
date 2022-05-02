//
// Created by Charles Du on 4/5/22.
//

#include "sTGC_3D_Formulation.h"
#include "geo_util.h"
#include <iostream>


sTGC_3D_Formulation::sTGC_3D_Formulation(const Matrix3Xd &rest_vertices, Matrix3Xd init_vertices, Matrix4Xi tets,
                                         const VectorXi &handles, const std::string &form, double alpha, double lambda1,
                                         double lambda2, double k, bool scale_rest_mesh,
                                         bool subtract_total_signed_content) :
        T(std::move(tets)), V(std::move(init_vertices)), subtract_total_signed_content(subtract_total_signed_content),
        skip_non_free_tets(false) {
    // compute freeI: indices of free vertices
    int nV = rest_vertices.cols();
    int vDim = 3;

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

    // experimental: compute free tet Indices
    // free tet: a tetrahedron with at least one free vertices
    int nT = T.cols();
    int n_free_tet = 0;
    std::vector<bool> free_tetQ(nT, false);
    for (int i = 0; i < nT; ++i) {
        int i1 = T(0,i);
        int i2 = T(1,i);
        int i3 = T(2,i);
        int i4 = T(3,i);
        if (freeQ[i1] || freeQ[i2] || freeQ[i3] || freeQ[i4]) {
            free_tetQ[i] = true;
            ++n_free_tet;
        }
    }
//    std::cout << "n_free_tet = " << n_free_tet << std::endl;
    free_tetI.resize(n_free_tet);
    ii = 0;
    for (int i = 0; i < nT; ++i) {
        if (free_tetQ[i]) {
            free_tetI[ii] = i;
            ++ii;
        }
    }


    // compute indexDict and F_free
    indexDict = VectorXi::Constant(nV, -1);
    for (auto i = 0; i < freeI.size(); ++i) {
        indexDict(freeI(i)) = i;
    }

    T_free.resize(T.rows(), T.cols());
    for (auto i = 0; i < T.cols(); ++i) {
        for (auto j = 0; j < T.rows(); ++j) {
            T_free(j, i) = indexDict(T(j, i));
        }
    }

    // (optional) scale rest mesh to have the same volume as init mesh
    double init_volume = compute_total_signed_mesh_volume(V,T);
    double rest_volume = compute_total_unsigned_volume(rest_vertices,T);
    double rest_scale = 1;
    if (scale_rest_mesh && handles.size() > 1 && init_volume > 0) {
        rest_scale = cbrt(init_volume / rest_volume);
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
    tgc.initialize(scaled_rest_vertices, T, form, alpha, lambda1, lambda2, k);
    if (skip_non_free_tets) {
        tgc.set_free_tetI(free_tetI);
    }

    // extract boundary triangles
    extract_mesh_boundary_triangles(T, boundary_triangles);

    // check if the boundary is fixed
    fixed_boundary = true;
    for (auto & boundary_triangle : boundary_triangles) {
        if (freeQ[boundary_triangle[0]] || freeQ[boundary_triangle[1]] || freeQ[boundary_triangle[2]])
        {
            fixed_boundary = false;
            break;
        }
    }
    boundary_signed_content = 0;
    if (fixed_boundary) {
        boundary_signed_content = init_volume;
    }
}

double sTGC_3D_Formulation::compute_energy(const VectorXd &x) {
    update_V(x);

    if (subtract_total_signed_content) {
        if (fixed_boundary) {
            return tgc.compute_total_generalized_content(V) - boundary_signed_content;
        } else {
            return tgc.compute_total_generalized_negative_content(V);
        }
    } else {
        return tgc.compute_total_generalized_content(V);
    }
}

void sTGC_3D_Formulation::update_V(const VectorXd &x) {
    int vDim = 3;
    for (auto i = 0; i < freeI.size(); ++i) {
        for (auto j = 0; j < vDim; ++j) {
            V(j, freeI(i)) = x[i * vDim + j];
        }
    }
}

double sTGC_3D_Formulation::compute_energy(const VectorXd &x, VectorXd &energy_list) {
    update_V(x);

    if (subtract_total_signed_content) {
        return tgc.compute_total_generalized_negative_content(V, energy_list);
    } else {
        return tgc.compute_total_generalized_content(V, energy_list);
    }
}

double sTGC_3D_Formulation::compute_energy_with_gradient(const VectorXd &x, VectorXd &grad) {
    update_V(x);

    // TGC
    Matrix3Xd tgc_grad;
    double tgc_energy = tgc.compute_total_generalized_content_with_gradient(V, tgc_grad);

    if (subtract_total_signed_content && !fixed_boundary) {
        //
        std::vector<Point3D> vertices(V.cols());
        for (int i = 0; i < V.cols(); ++i) {
            vertices[i] = V.col(i);
        }
        Matrix3Xd total_signed_content_grad;
        total_signed_content_grad.setZero(3, vertices.size());
        double total_signed_content = compute_total_signed_volume_with_gradient(
                vertices, boundary_triangles, total_signed_content_grad);

        // full gradient
        Matrix3Xd m_grad = tgc_grad - total_signed_content_grad;

        // gradient of free vertices
        grad.resize(freeI.size() * m_grad.rows());
        for (int i = 0; i < freeI.size(); ++i) {
            grad(3*i) = m_grad(0,freeI(i));
            grad(3*i+1) = m_grad(1,freeI(i));
            grad(3*i+2) = m_grad(2,freeI(i));
        }

        // energy
        return tgc_energy - total_signed_content;
    }
    else {
        // gradient of free vertices
        grad.resize(freeI.size() * tgc_grad.rows());
        for (int i = 0; i < freeI.size(); ++i) {
            grad(3*i) = tgc_grad(0,freeI(i));
            grad(3*i+1) = tgc_grad(1,freeI(i));
            grad(3*i+2) = tgc_grad(2,freeI(i));
        }

        // energy
        if (subtract_total_signed_content) { // the boundary must be fixed
            return tgc_energy - boundary_signed_content;
        } else {
            return tgc_energy;
        }
    }
}

double sTGC_3D_Formulation::compute_energy_with_gradient_projectedHessian(const VectorXd &x, VectorXd &energy_list,
                                                                          VectorXd &grad, SpMat &Hess) {
    update_V(x);

    // TGC energy, TGC full gradient,
    // and PSD projected Hessian of (TGC - total signed content) on free vertices
    VectorXd generalized_content_list;
    Matrix3Xd tgc_grad;
    double tgc_energy = tgc.compute_total_generalized_content_with_gradient_and_sTGC_projectedHessian(V, freeI,
                                                                                                      T_free,
                                                                                                      generalized_content_list,
                                                                                                      tgc_grad,
                                                                                                      Hess);

    if (!subtract_total_signed_content) {
        energy_list = generalized_content_list;
    } else {
        // signed volumes
        VectorXd signed_volume_list;
        compute_signed_tet_volumes(V, T, signed_volume_list);

        // fill the energy decomposition into energy_list
//        energy_list.resize(generalized_content_list.size());
//        for (int i = 0; i < generalized_content_list.size(); ++i) {
//            energy_list(i) = generalized_content_list(i) - signed_volume_list(i);
//        }
        energy_list = generalized_content_list - signed_volume_list;
    }

    if (subtract_total_signed_content && !fixed_boundary) {
        //
        std::vector<Point3D> vertices(V.cols());
        for (int i = 0; i < V.cols(); ++i) {
            vertices[i] = V.col(i);
        }
        Matrix3Xd total_signed_content_grad;
        total_signed_content_grad.setZero(3, vertices.size());
        double total_signed_content = compute_total_signed_volume_with_gradient(
                vertices, boundary_triangles, total_signed_content_grad);

        // full gradient
        Matrix3Xd m_grad = tgc_grad - total_signed_content_grad;

        // gradient of free vertices
        grad.resize(freeI.size() * m_grad.rows());
        for (int i = 0; i < freeI.size(); ++i) {
            grad(3*i) = m_grad(0,freeI(i));
            grad(3*i+1) = m_grad(1,freeI(i));
            grad(3*i+2) = m_grad(2,freeI(i));
        }

        // energy
        return tgc_energy - total_signed_content;
    }
    else {
        // gradient of free vertices
        grad.resize(freeI.size() * tgc_grad.rows());
        for (int i = 0; i < freeI.size(); ++i) {
            grad(3*i) = tgc_grad(0,freeI(i));
            grad(3*i+1) = tgc_grad(1,freeI(i));
            grad(3*i+2) = tgc_grad(2,freeI(i));
        }

        // energy
        if (subtract_total_signed_content) { // the boundary must be fixed
            return tgc_energy - boundary_signed_content;
        } else {
            return tgc_energy;
        }
    }
}

