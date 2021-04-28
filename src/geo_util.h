//
// Created by Charles Du on 4/21/21.
//

#ifndef TLC_GEO_UTIL_H
#define TLC_GEO_UTIL_H

#include <Eigen/Core>
#include <vector>

#include <autodiff/reverse.hpp>
#include <autodiff/reverse/eigen.hpp>

typedef Eigen::Matrix<autodiff::var,2,1> Point;

typedef Eigen::Matrix<autodiff::var,2,1> VarVec2d;

typedef autodiff::var var;

// return the vector after rotating 90 degree counter clockwise
VarVec2d rotate_90deg(const VarVec2d &vec);

// compute the angle mod 2pi
// return: angle in [0, 2pi)
var angle_mod_2PI(var a);

// compute angle from positive x-axis to vec
// return: angle in [0, 2pi)
var compute_vector_angle(const Point &p);

// compute the counter clockwise rotation angle from a1 to a2
// return: angle in [0, 2pi)
var compute_rotation_angle(var a1, var a2);

// for a ray sweeping counter-clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
bool is_angle_between_ccw(var a, var a1, var a2);

// for a ray sweeping clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
bool is_angle_between_cw(var a, var a1, var a2);

// test if angle a is between a1(included) and a2(not included)
bool is_angle_between(var a, var a1, var a2);

// compute total signed area of a polygon
var compute_total_signed_area(const std::vector<Point> &vertices,
                                 const std::vector<std::pair<size_t,size_t>> &edges);


#endif //TLC_GEO_UTIL_H

