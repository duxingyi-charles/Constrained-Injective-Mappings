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
    // parent edge index
    size_t parent_id;
    // circular arc geometry
    Circular_Arc arc;
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

    // decompose 2D plane into arrangement cells
    // input: subdivided polyArc {pts, pEdges}
    // output: cells, each cell represented by a list of boundary half-edges
    // notes:
    // * index convention for half-edge:
    // ** i : index of original edge e
    // ** 2i+1: same direction as e, on the left side of e
    // ** 2i  : opposite direction of e, on the right side of e
    static void decompose_into_cells(const std::vector<Point> &pts, const std::vector<SubArc_Edge> &pEdges,
                                     // output
                                     std::vector<std::vector<size_t>> &cells);

private:
    // trace half-edges into chains
    // input: next_hEdge[i] is the index of the half-edge after half-edge i.
    // output: chains
    // - each chain is a connected component of half-edges
    static void trace_chains(const std::vector<size_t> &next_hEdge,
                      // output
                      std::vector<std::vector<size_t>> &chains);

};


#endif //TLC_ARRANGEMENT_H
