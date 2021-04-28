//
// Created by Charles Du on 4/21/21.
//

#include "geo_util.h"

Eigen::Vector2d rotate_90deg(const Eigen::Vector2d &vec) {
//    return Eigen::Vector2d(-(vec.y()), vec.x());
    return {-(vec.y()), vec.x()};
}


double angle_mod_2PI(double a) {
    if (a >= 0) {
        return fmod(a, 2*M_PI);
    } else {
        return 2*M_PI + fmod(a, 2*M_PI);
    }
}

double compute_vector_angle(const Eigen::Vector2d &p) {
    double x = p.x(), y = p.y();
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

double compute_rotation_angle(double a1, double a2) {
    double angle1 = angle_mod_2PI(a1);
    double angle2 = angle_mod_2PI(a2);
    if (angle1 <= angle2) {
        return angle2 - angle1;
    } else {
        return angle2 + 2*M_PI - angle1;
    }
}

bool is_angle_between_ccw(double a, double a1, double a2) {
    double angle = angle_mod_2PI(a);
    double angle1 = angle_mod_2PI(a1);
    double angle2 = angle_mod_2PI(a2);

    if (angle1 <= angle2) {
        return (angle >= angle1) && (angle < angle2);
    } else {
        return (angle < angle2) || (angle >= angle1);
    }
}

bool is_angle_between_cw(double a, double a1, double a2) {
    double angle = angle_mod_2PI(a);
    double angle1 = angle_mod_2PI(a1);
    double angle2 = angle_mod_2PI(a2);

    if (angle2 <= angle1) {
        return (angle2 < angle) && (angle <= angle1);
    } else {
        return (angle <= angle1) || (angle > angle2);
    }
}

bool is_angle_between(double a, double a1, double a2) {
    if (a1 < a2) {
        return is_angle_between_ccw(a,a1,a2);
    } else {
        return is_angle_between_cw(a,a1,a2);
    }
}

double compute_total_signed_area(const std::vector<Point> &vertices,
                                 const std::vector<std::pair<size_t, size_t>> &edges) {
    double area = 0;
    for (const auto & e : edges) {
        auto i = e.first;
        auto j = e.second;
        area += vertices[i].x() * vertices[j].y() - vertices[i].y() * vertices[j].x();
    }
    area /= 2;

    return area;
}

