//
// Created by Charles Du on 4/28/21.
//

#include "Arc_Occupancy.h"
#include "Arrangement.h"

// explicit instantiations
template class Arc_Occupancy<double>;
//

template <typename Scalar>
Scalar Arc_Occupancy<Scalar>::compute_arc_occupancy(const std::vector<Point<Scalar>> &vertices,
                                            const std::vector<std::pair<size_t, size_t>> &edges) const {

    return compute_arc_occupancy(vertices, edges, param_theta);
}

template <typename Scalar>
Scalar Arc_Occupancy<Scalar>::compute_arc_loop_area(const std::vector<Point<Scalar>> &pts,
                                                    const std::vector<SubArc_Edge<Scalar>> &edges)
{
    std::vector<std::pair<size_t, size_t>> raw_edges;
    raw_edges.reserve(edges.size());
    for (const auto & e : edges) {
        raw_edges.emplace_back(e.id1, e.id2);
    }

    // signed area of the polygon
    auto polygon_area = compute_total_signed_area(pts, raw_edges);

    // area of arc segments
    Scalar arc_seg_area = 0;
    for (const auto & e: edges) {
        arc_seg_area += e.arc.get_segment_area();
    }

    //
    return polygon_area + arc_seg_area;
}

template <typename Scalar>
Scalar Arc_Occupancy<Scalar>::compute_arc_occupancy(const std::vector<Point<Scalar>> &vertices,
                                            const std::vector<std::pair<size_t, size_t>> &edges, Scalar theta) {

    // create arc edges
    std::vector<Arc_Edge<Scalar>> arc_edges;
    arc_edges.reserve(edges.size());

    for (const auto & e : edges) {
        arc_edges.emplace_back(Arc_Edge{e.first, e.second,
                                        Circular_Arc(vertices[e.first], vertices[e.second], theta)});
    }

    // compute arrangement and winding numbers
    std::vector<Point<Scalar>> pts;
    std::vector<SubArc_Edge<Scalar>> pEdges;
    std::vector<bool> is_intersection_point;
    std::vector<std::vector<SubArc_Edge<Scalar>>> edges_of_cell;
    std::vector<int> windings;
    Arrangement<Scalar>::compute_arrangement(vertices,arc_edges,
                                     pts,pEdges,is_intersection_point,edges_of_cell,windings);

    // compute occupancy
    Scalar occupancy = 0;
    for (int i = 0; i < windings.size(); ++i) {
        if (windings[i] > 0) {
            occupancy += compute_arc_loop_area(pts, edges_of_cell[i]);
        }
    }

    //
    return occupancy;
}
