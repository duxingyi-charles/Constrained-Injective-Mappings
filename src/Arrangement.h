//
// Created by Charles Du on 4/26/21.
//

#ifndef TLC_ARRANGEMENT_H
#define TLC_ARRANGEMENT_H

#include "Circular_Arc.h"

struct Arc_Edge {
    // first vertex index
    size_t id1;
    // second vertex index
    size_t id2;
    // circular arc geometry
    Circular_Arc arc;
};

struct SubArc_Edge {
    // first vertex index
    size_t id1;
    // second vertex index
    size_t id2;
    // circular arc angle
    double arc_angle;
    // parent edge index
    size_t parent_id;
};

class Arrangement {
public:
    Arrangement() = default;
    ~Arrangement() = default;

    // compute self-arrangement of a circular arc loop
    // save result as internal states
    void compute_arrangement(const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges);

    // subdivide input arcs by their intersections
    static void subdivide_polyArc_by_intersection(const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges,
                                             // output
                                             std::vector<Point> &pts,
                                             std::vector<SubArc_Edge> &pEdges,
                                             std::vector<bool> &is_intersection_point);

private:


};


#endif //TLC_ARRANGEMENT_H
