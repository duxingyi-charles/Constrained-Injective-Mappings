//
// Created by Charles Du on 4/21/21.
//

#include "geo_util.h"


double get_vector_angle(const Eigen::Vector2d &vec) {
    double x = vec.x(), y = vec.y();
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


bool is_angle_between_ccw(double a, double a1, double a2) {
    double angle = fmod(a, 2*M_PI);
    double angle1 = fmod(a1, 2*M_PI);
    double angle2 = fmod(a2, 2*M_PI);

    if (angle1 <= angle2) {
        return (angle >= angle1) && (angle < angle2);
    } else {
        return (angle < angle2) || (angle >= angle1);
    }
}

bool is_angle_between_cw(double a, double a1, double a2) {
    double angle = fmod(a, 2*M_PI);
    double angle1 = fmod(a1, 2*M_PI);
    double angle2 = fmod(a2, 2*M_PI);

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

