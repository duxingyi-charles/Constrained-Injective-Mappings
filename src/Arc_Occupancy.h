//
// Created by Charles Du on 4/28/21.
//

#ifndef TLC_ARC_OCCUPANCY_H
#define TLC_ARC_OCCUPANCY_H

#include "Arrangement.h"

template <typename Scalar>
class Arc_Occupancy {
public:
    Arc_Occupancy() = default;
    explicit Arc_Occupancy(Scalar t) : param_theta(t) {};
    ~Arc_Occupancy() = default;

    // compute arc occupancy (static)
    static Scalar compute_arc_occupancy(const std::vector<Point<Scalar>> &vertices,
                                        const std::vector<std::pair<size_t,size_t>> &edges,
                                        Scalar theta);

    // compute arc occupancy
    Scalar compute_arc_occupancy(const std::vector<Point<Scalar>> &vertices,
                                 const std::vector<std::pair<size_t,size_t>> &edges) const;

    // compute signed area of an arc loop
    static Scalar compute_arc_loop_area(const std::vector<Point<Scalar>> &pts,
                                        const std::vector<SubArc_Edge<Scalar>> &edges);

private:
    // arc angle
    Scalar param_theta;

};


#endif //TLC_ARC_OCCUPANCY_H
