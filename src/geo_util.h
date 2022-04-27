//
// Created by Charles Du on 4/21/21.
//

#ifndef TLC_GEO_UTIL_H
#define TLC_GEO_UTIL_H

#include <Eigen/Core>
#include <vector>
#include <array>

#include "corecrt_math_defines.h"

typedef Eigen::Vector2d Point;
typedef Eigen::Vector3d Point3D;

// return the vector after rotating 90 degree counter clockwise
Eigen::Vector2d rotate_90deg(const Eigen::Vector2d &vec);

// compute the angle mod 2pi
// return: angle in [0, 2pi)
double angle_mod_2PI(double a);

// compute angle from positive x-axis to vec
// return: angle in [0, 2pi)
double compute_vector_angle(const Point &p);

// compute the counter clockwise rotation angle from a1 to a2
// return: angle in [0, 2pi)
double compute_rotation_angle(double a1, double a2);

// for a ray sweeping counter-clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
bool is_angle_between_ccw(double a, double a1, double a2);

// for a ray sweeping clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
bool is_angle_between_cw(double a, double a1, double a2);

// test if angle a is between a1(included) and a2(not included)
bool is_angle_between(double a, double a1, double a2);

// compute total signed area of a polygon
double compute_total_signed_area(const std::vector<Point> &vertices,
                                 const std::vector<std::pair<size_t,size_t>> &edges);

// compute total signed volume of a closed surface
double compute_total_signed_volume(const std::vector<Point3D> &vertices,
                                   const std::vector<std::array<int,3>> &triangles);

// compute total signed area of a polygon and its gradient wrt. polygon vertices
double compute_total_signed_area_with_gradient(const std::vector<Point> &vertices,
                                               const std::vector<std::pair<size_t,size_t>> &edges,
                                               Eigen::Matrix2Xd &dArea_dv);

// compute total signed volume of a closed surface and its gradient wrt. surface vertices
double compute_total_signed_volume_with_gradient(const std::vector<Point3D> &vertices,
                                                 const std::vector<std::array<int,3>> &triangles,
                                                 Eigen::Matrix3Xd &dVolume_dv);

// compute triangle area using (robust) Heron's formula
// input: squared edge lengths of the triangle
double compute_Heron_tri_area(double d1, double d2, double d3);

// compute triangle aspect ratio
// input: squared edge lengths of the triangle
double compute_tri_aspect_ratio(double d1, double d2, double d3);

// computed total signed area of a triangle mesh
double compute_total_signed_mesh_area(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F);

// compute total signed volume of a tetrahedron mesh
double compute_total_signed_mesh_volume(const Eigen::Matrix3Xd &V,const Eigen::Matrix4Xi &Tets);

// compute minimum signed area of a triangle mesh
double compute_min_signed_mesh_area(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F);

// compute minimum signed volume of a tetrahedron mesh
double compute_min_signed_mesh_volume(const Eigen::Matrix3Xd &V, const Eigen::Matrix4Xi &Tets);

// compute signed area of all triangles in a mesh
void compute_signed_tri_areas(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F,
                                Eigen::VectorXd &areaList);

// compute signed volume of all tetrahedrons in a mesh
void compute_signed_tet_volumes(const Eigen::Matrix3Xd &V, const Eigen::Matrix4Xi &Tets,
                              Eigen::VectorXd &volumeList);

// compute signed area of a triangle
// input: 2D coordinates of 3 points of triangle
// return: signed area of the triangle
double compute_tri_signed_area(const Point &p1, const Point &p2, const Point &p3);

// compute signed volume of a tetrahedron
// input: 3D coordinates of 4 points of tetrahedron
// return: signed volume of the tetrahedron
double compute_tet_signed_volume(const Point3D &p1, const Point3D &p2, const Point3D &p3, const Point3D &p4);

// compute signed volume of a tetrahedron and its gradient w.r.t. vertices
// input: 3D coordinates of 4 points of tetrahedron
// return: signed volume of the tetrahedron, and its gradient
double compute_tet_signed_volume_with_gradient(const Point3D &p1, const Point3D &p2, const Point3D &p3, const Point3D &p4,
                                               Eigen::Matrix3Xd &grad);


// compute total unsigned area of a triangle mesh
double compute_total_unsigned_area(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F);

// compute total unsigned volume of a tetrahedron mesh
double compute_total_unsigned_volume(const Eigen::Matrix3Xd &V, const Eigen::Matrix4Xi &Tets);

//
void compute_squared_edge_Length(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
                                 Eigen::Matrix3Xd &D);

//
void extract_mesh_boundary_edges(const Eigen::Matrix3Xi &faces,
                                 std::vector<std::pair<size_t,size_t>> &boundary_edges);

//
void extract_mesh_boundary_triangles(const Eigen::Matrix4Xi &tets,
                                 std::vector<std::array<int,3>>& boundary_triangles);

// compute angle between from vector1 to vector2
// return: angle in (-pi, pi]
double compute_vec_vec_angle(const Eigen::Vector2d &vec1, const Eigen::Vector2d &vec2);



// compute indices of over-winded interior vertices
// input: 2D triangle mesh (V, F); is_boundary_vertex
// output: indices of winded interior vertices
void compute_winded_interior_vertices(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F,
                                      const std::vector<bool> &is_boundary_vertex,
                                      std::vector<size_t> &winded_vertices);

#endif //TLC_GEO_UTIL_H

