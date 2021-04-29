//
// Created by Charles Du on 4/28/21.
//

#ifndef TLC_ARC_OCCUPANCY_H
#define TLC_ARC_OCCUPANCY_H

#include "Arrangement.h"

class Arc_Occupancy {
public:
    explicit Arc_Occupancy(var t) : param_theta(t) {};
    ~Arc_Occupancy() = default;

    // compute arc occupancy (static)
    static var compute_arc_occupancy(const std::vector<Point> &vertices,
                                        const std::vector<std::pair<size_t,size_t>> &edges,
                                        var theta);

    // compute arc occupancy
    var compute_arc_occupancy(const std::vector<Point> &vertices,
                                 const std::vector<std::pair<size_t,size_t>> &edges) const;

    // compute signed area of an arc loop
    static var compute_arc_loop_area(const std::vector<Point> &pts, const std::vector<SubArc_Edge> &edges);

private:
    // arc angle
    var param_theta;

};


#endif //TLC_ARC_OCCUPANCY_H
