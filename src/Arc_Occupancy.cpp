//
// Created by Charles Du on 4/28/21.
//

#include "Arc_Occupancy.h"
#include "Arrangement.h"

var Arc_Occupancy::compute_arc_occupancy(const std::vector<Point> &vertices,
                                            const std::vector<std::pair<size_t, size_t>> &edges) const {

    return compute_arc_occupancy(vertices, edges, param_theta);
}

var Arc_Occupancy::compute_arc_loop_area(const std::vector<Point> &pts, const std::vector<SubArc_Edge> &edges)
{
    std::vector<std::pair<size_t, size_t>> raw_edges;
    raw_edges.reserve(edges.size());
    for (const auto & e : edges) {
        raw_edges.emplace_back(e.id1, e.id2);
    }

    // signed area of the polygon
    var polygon_area = compute_total_signed_area(pts, raw_edges);

    // area of arc segments
    var arc_seg_area = 0;
    for (const auto & e: edges) {
        arc_seg_area += e.arc.get_segment_area();
    }

    //
    return polygon_area + arc_seg_area;
}

var Arc_Occupancy::compute_arc_occupancy(const std::vector<Point> &vertices,
                                            const std::vector<std::pair<size_t, size_t>> &edges, var theta) {

    // create arc edges
    std::vector<Arc_Edge> arc_edges;
    arc_edges.reserve(edges.size());

    for (const auto & e : edges) {
        arc_edges.emplace_back(Arc_Edge{e.first, e.second,
                                        Circular_Arc(vertices[e.first], vertices[e.second], theta)});
    }

    // compute arrangement and winding numbers
    std::vector<Point> pts;
    std::vector<SubArc_Edge> pEdges;
    std::vector<bool> is_intersection_point;
    std::vector<std::vector<SubArc_Edge>> edges_of_cell;
    std::vector<int> windings;
    Arrangement::compute_arrangement(vertices,arc_edges,
                                     pts,pEdges,is_intersection_point,edges_of_cell,windings);

    // compute occupancy
    var occupancy = 0;
    for (int i = 0; i < windings.size(); ++i) {
        if (windings[i] > 0) {
            occupancy += compute_arc_loop_area(pts, edges_of_cell[i]);
        }
    }

    //
    return occupancy;
}

var Arc_Occupancy::compute_arc_occupancy(const Eigen::VectorXvar &vertices,
                          const std::vector<std::pair<size_t,size_t>> &edges) const
{
    // convert vertices to a vector of Point
    std::vector<Point> verts;
    auto n_vert = vertices.size() /2;
    verts.reserve(n_vert);

    for (int i = 0; i < n_vert; ++i) {
        verts.emplace_back(vertices[2*i], vertices[2*i+1]);
    }

    //
    return compute_arc_occupancy(verts, edges, param_theta);
}
