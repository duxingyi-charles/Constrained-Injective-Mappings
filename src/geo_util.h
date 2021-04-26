//
// Created by Charles Du on 4/21/21.
//

#ifndef TLC_GEO_UTIL_H
#define TLC_GEO_UTIL_H

#include <Eigen/Core>

typedef Eigen::Vector2d Point;

// compute angle from positive x-axis to vec
// return: angle in [0, 2pi)
double get_vector_angle(const Point &p);

// for a ray sweeping counter-clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
bool is_angle_between_ccw(double a, double a1, double a2);

// for a ray sweeping clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
bool is_angle_between_cw(double a, double a1, double a2);

// test if angle a is between a1(included) and a2(not included)
bool is_angle_between(double a, double a1, double a2);

#endif //TLC_GEO_UTIL_H

