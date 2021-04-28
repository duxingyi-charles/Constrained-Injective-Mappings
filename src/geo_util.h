//
// Created by Charles Du on 4/21/21.
//

#ifndef TLC_GEO_UTIL_H
#define TLC_GEO_UTIL_H

#include <Eigen/Core>
#include <vector>

#include <autodiff/reverse.hpp>
#include <autodiff/reverse/eigen.hpp>

//typedef Eigen::Vector2d Point;
//typedef Eigen::Matrix<autodiff::var, 2, 1> Point;
template <typename Scalar>
using Point = Eigen::Matrix<Scalar,2,1>;

template <typename Scalar>
using Vec2D = Eigen::Matrix<Scalar,2,1>;

// return the vector after rotating 90 degree counter clockwise
template <typename Scalar>
Vec2D<Scalar> rotate_90deg(const Vec2D<Scalar> &vec);

// compute the angle mod 2pi
// return: angle in [0, 2pi)
template <typename Scalar>
Scalar angle_mod_2PI(Scalar a);

// compute angle from positive x-axis to vec
// return: angle in [0, 2pi)
template <typename Scalar>
Scalar compute_vector_angle(const Point<Scalar> &p);

// compute the counter clockwise rotation angle from a1 to a2
// return: angle in [0, 2pi)
template <typename Scalar>
Scalar compute_rotation_angle(Scalar a1, Scalar a2);

// for a ray sweeping counter-clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
template <typename Scalar>
bool is_angle_between_ccw(Scalar a, Scalar a1, Scalar a2);

// for a ray sweeping clockwise from angle a1(included) to a2(not included),
// test if the ray will encounter angle a
template <typename Scalar>
bool is_angle_between_cw(Scalar a, Scalar a1, Scalar a2);

// test if angle a is between a1(included) and a2(not included)
template <typename Scalar>
bool is_angle_between(Scalar a, Scalar a1, Scalar a2);

// compute total signed area of a polygon
//double
template<typename Scalar>
Scalar compute_total_signed_area(const std::vector<Point<Scalar>> &vertices,
                                 const std::vector<std::pair<size_t,size_t>> &edges);

// implementation

template <typename Scalar>
Vec2D<Scalar> rotate_90deg(const Vec2D<Scalar> &vec) {
//    return Eigen::Vector2d(-(vec.y()), vec.x());
    return {-(vec.y()), vec.x()};
}

template <typename Scalar>
Scalar angle_mod_2PI(Scalar a) {
    if (a >= 0) {
        return fmod(a, 2*M_PI);
    } else {
        return 2*M_PI + fmod(a, 2*M_PI);
    }
}

template <typename Scalar>
Scalar compute_vector_angle(const Point<Scalar> &p) {
    Scalar x = p.x(), y = p.y();
    if (x == 0) {
        if (y >=0) {
            return M_PI/2;
        } else {
            return 3*M_PI/2;
        }
    } else {
        if (x > 0) {
            if (y >= 0) {
                return asin(y/sqrt(x*x+y*y));
            } else {
                return 2*M_PI + asin(y/sqrt(x*x+y*y));
            }
        } else { // x < 0
            return M_PI - asin(y/sqrt(x*x+y*y));
        }
    }
}

template <typename Scalar>
Scalar compute_rotation_angle(Scalar a1, Scalar a2) {
    Scalar angle1 = angle_mod_2PI(a1);
    Scalar angle2 = angle_mod_2PI(a2);
    if (angle1 <= angle2) {
        return angle2 - angle1;
    } else {
        return angle2 + 2*M_PI - angle1;
    }
}

template <typename Scalar>
bool is_angle_between_ccw(Scalar a, Scalar a1, Scalar a2) {
    Scalar angle = angle_mod_2PI(a);
    Scalar angle1 = angle_mod_2PI(a1);
    Scalar angle2 = angle_mod_2PI(a2);

    if (angle1 <= angle2) {
        return (angle >= angle1) && (angle < angle2);
    } else {
        return (angle < angle2) || (angle >= angle1);
    }
}

template <typename Scalar>
bool is_angle_between_cw(Scalar a, Scalar a1, Scalar a2) {
    Scalar angle = angle_mod_2PI(a);
    Scalar angle1 = angle_mod_2PI(a1);
    Scalar angle2 = angle_mod_2PI(a2);

    if (angle2 <= angle1) {
        return (angle2 < angle) && (angle <= angle1);
    } else {
        return (angle <= angle1) || (angle > angle2);
    }
}

template <typename Scalar>
bool is_angle_between(Scalar a, Scalar a1, Scalar a2) {
    if (a1 < a2) {
        return is_angle_between_ccw(a,a1,a2);
    } else {
        return is_angle_between_cw(a,a1,a2);
    }
}

template<typename Scalar>
Scalar
compute_total_signed_area(const std::vector<Point<Scalar>> &vertices, const std::vector<std::pair<size_t, size_t>> &edges)
{
    Scalar area = 0;
    for (const auto & e : edges) {
        auto i = e.first;
        auto j = e.second;
        area += vertices[i].x() * vertices[j].y() - vertices[i].y() * vertices[j].x();
    }
    area /= 2;

    return area;
}

#endif //TLC_GEO_UTIL_H

