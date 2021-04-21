//
// Created by Charles Du on 4/21/21.
//

#include "geo_util.h"

// compute angle from positive x-axis to vec
// return: angle in [0, 2pi)
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
