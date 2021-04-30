//
// Created by Charles Du on 4/28/21.
//

#ifndef TLC_ARC_OCCUPANCY_H
#define TLC_ARC_OCCUPANCY_H

#include "Arrangement.h"

class Arc_Occupancy {
public:
    explicit Arc_Occupancy(double t) : param_theta(t) {};
    ~Arc_Occupancy() = default;

    // compute arc occupancy (static)
    static double compute_arc_occupancy(const std::vector<Point> &vertices,
                                        const std::vector<std::pair<size_t,size_t>> &edges,
                                        double theta);

    // compute arc occupancy and its gradient wrt. input vertices
    static double compute_arc_occupancy_with_gradient(const std::vector<Point> &vertices,
                                                      const std::vector<std::pair<size_t,size_t>> &edges,
                                                      double theta,
                                                      // output
                                                      Eigen::Matrix2Xd &grad);

    // compute arc occupancy and its gradient wrt. input vertices
    double compute_arc_occupancy_with_gradient(const std::vector<Point> &vertices,
                                               const std::vector<std::pair<size_t,size_t>> &edges,
                                               // output
                                               Eigen::Matrix2Xd &grad);

    // compute arc occupancy
    double compute_arc_occupancy(const std::vector<Point> &vertices,
                                 const std::vector<std::pair<size_t,size_t>> &edges) const;

    // compute signed area of an arc loop
    static double compute_arc_loop_area(const std::vector<Point> &pts, const std::vector<SubArc_Edge> &edges);

    // compute arc curve segment area
    double compute_arc_curve_segment_area(const std::vector<Point> &vertices,
                                 const std::vector<std::pair<size_t,size_t>> &edges) const;

    double compute_arc_curve_segment_area_with_gradient(const std::vector<Point> &vertices,
                                                        const std::vector<std::pair<size_t,size_t>> &edges,
                                                        //output
                                                        Eigen::Matrix2Xd &grad) const;


private:
    // arc angle
    double param_theta;

    // compute signed area of an arc loop and update the gradient
    static double compute_arc_loop_area_with_gradient(const std::vector<Point> &pts, const std::vector<SubArc_Edge> &edges,
                                                      Eigen::Matrix2Xd &dArea_dp, Eigen::VectorXd &dArea_dr2);

};


#endif //TLC_ARC_OCCUPANCY_H
