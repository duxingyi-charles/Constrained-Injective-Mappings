//
// Created by Charles Du on 4/29/21.
//

#include "Total_Lifted_Content.h"

#include <utility>

void Total_Lifted_Content::initialize(const Eigen::MatrixXd &rest_vertices, Eigen::Matrix3Xi faces,
                                           const std::string &form, double alpha)
{
    //
    F = faces;

    // compute restD
    if (form == "harmonic") {
        compute_squared_edge_Length(rest_vertices, F, restD);
        restD *= alpha;
    } else // Tutte form
    {
        restD = Eigen::MatrixXd::Constant(3, F.cols(), alpha);
    }
}



double Total_Lifted_Content::compute_total_lifted_content(const Matrix2Xd &vertices) const {

    VectorXd energyList(F.cols());
    int vDim = 2;
    int simplex_size = 3; //triangle

    for (auto i = 0; i < F.cols(); ++i) {
        Matrix2Xd vert(vDim, simplex_size);
        vert.col(0) = vertices.col(F(0, i));
        vert.col(1) = vertices.col(F(1, i));
        vert.col(2) = vertices.col(F(2, i));

        const Vector3d &r = restD.col(i);

        energyList(i) = compute_lifted_TriArea(vert, r);
    }

    return energyList.sum();
}

double Total_Lifted_Content::compute_lifted_TriArea(const Matrix2Xd &vert, const Vector3d &r) {
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto e1 = v2 - v3;
    auto e2 = v3 - v1;
    auto e3 = v1 - v2;
    double d1 = e1.squaredNorm() + r(0);
    double d2 = e2.squaredNorm() + r(1);
    double d3 = e3.squaredNorm() + r(2);
    return compute_Heron_tri_area(d1, d2, d3);
}

double
Total_Lifted_Content::compute_total_lifted_content_with_gradient(const Matrix2Xd &vertices, Matrix2Xd &grad) const
{
    int vDim = 2;
    double energy = 0.0;
    grad = Matrix2Xd::Zero(2, vertices.cols());

    for (auto i = 0; i < F.cols(); ++i) {
        int i1, i2, i3;
        i1 = F(0, i);
        i2 = F(1, i);
        i3 = F(2, i);

        MatrixXd vert(vDim, 3);
        vert.col(0) = vertices.col(i1);
        vert.col(1) = vertices.col(i2);
        vert.col(2) = vertices.col(i3);
        Vector3d r = restD.col(i);


        Matrix2Xd g;
        energy += compute_lifted_TriArea_with_gradient(vert,r,g);

        grad.col(i1) += g.col(0);
        grad.col(i2) += g.col(1);
        grad.col(i3) += g.col(2);
    }

    return energy;
}

double
Total_Lifted_Content::compute_lifted_TriArea_with_gradient(const Matrix2Xd &vert, const Vector3d &r,
                                                           Matrix2Xd &grad)
{
    auto v1 = vert.col(0);
    auto v2 = vert.col(1);
    auto v3 = vert.col(2);
    auto e1 = v2 - v3;
    auto e2 = v3 - v1;
    auto e3 = v1 - v2;
    double d1 = e1.squaredNorm() + r(0);
    double d2 = e2.squaredNorm() + r(1);
    double d3 = e3.squaredNorm() + r(2);

    //
    double area = compute_Heron_tri_area(d1, d2, d3);

    //
    double g1 = d2 + d3 - d1;
    double g2 = d3 + d1 - d2;
    double g3 = d1 + d2 - d3;

    //
    auto ge1 = g1 * e1;
    auto ge2 = g2 * e2;
    auto ge3 = g3 * e3;

    //note: grad has the same dimension as vert
    grad.resize(vert.rows(), vert.cols());
    grad.col(0) = ge3 - ge2;
    grad.col(1) = ge1 - ge3;
    grad.col(2) = ge2 - ge1;

    double s = 1 / (8 * area);
    grad *= s;

    return area;
}


