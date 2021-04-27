//
// Created by Charles Du on 4/26/21.
//

#include <cassert>

#include <algorithm>
#include "Arrangement.h"

void Arrangement::compute_arrangement(const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges)
{
    // step 1: subdivide input arcs by their intersections
    std::vector<Point> pts;
    std::vector<SubArc_Edge> pEdges;
    std::vector<bool> is_intersection_point;
    subdivide_polyArc_by_intersection(vertices, edges, pts, pEdges,is_intersection_point);

    // step 2: get all arrangement cells
    std::vector<std::vector<size_t>> cells;
    decompose_into_cells(pts, pEdges, cells);

    // todo: step 3: compute winding numbers for each cell

    // save results


}


void Arrangement::subdivide_polyArc_by_intersection(
        const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges,
        std::vector<Point> &pts, std::vector<SubArc_Edge> &pEdges,
        std::vector<bool> &is_intersection_point)
{
    // init
    pts = vertices;
    pEdges.clear();

    is_intersection_point.clear();
    is_intersection_point.resize(pts.size(), false);

    typedef std::pair<size_t, double> PID_Angle_Pair;
    std::vector<std::vector<PID_Angle_Pair>> edge_intersection_list(edges.size()); // record intersections on each input arc edge

    // filter potentially intersecting arcs by bounding box
    std::vector<Rectangle> bbox_list;
    bbox_list.resize(edges.size());
    for (const auto& e : edges) {
        bbox_list.emplace_back(e.arc.get_bounding_box());
    }

    std::vector<std::pair<int,int>> potential_pair_list;
    for (int i = 0; i < edges.size(); ++i) {
        for (int j = i+1; j < edges.size(); ++j) {
            if (Rectangle::is_intersect(bbox_list[i], bbox_list[j])) {
                potential_pair_list.emplace_back(i, j);
            }
        }
    }

    // compute intersections
    for (const auto & pair : potential_pair_list) {
        int i = pair.first;
        int j = pair.second;
        std::vector<Intersection_Point> result;
        Circular_Arc::compute_intersection(edges[i].arc, edges[j].arc, result);
        if (!result.empty()) {
            for (const auto & intersection : result) {
                pts.emplace_back(intersection.location);
                is_intersection_point.push_back(true);
                edge_intersection_list[i].emplace_back(pts.size()-1, intersection.angle1);
                edge_intersection_list[j].emplace_back(pts.size()-1, intersection.angle2);
                // todo: intersectInfo

            }
        }
    }

    // subdivide input arcs
    auto PID_angle_less_than = [](const PID_Angle_Pair &left, const PID_Angle_Pair &right)
    { return left.second < right.second; };
    auto PID_angle_greater_than = [](const PID_Angle_Pair &left, const PID_Angle_Pair &right)
    { return left.second > right.second; };

    for (size_t ei = 0; ei < edges.size(); ++ei) {
        if (edge_intersection_list[ei].empty()) {
            // no intersection point in ei, copy the input edge
            pEdges.emplace_back(SubArc_Edge{edges[ei].id1, edges[ei].id2, ei, edges[ei].arc});
        } else {
            double theta1 = edges[ei].arc.get_start_angle();
            double theta2 = edges[ei].arc.get_end_angle();
            double theta  = edges[ei].arc.get_arc_angle();

            std::vector<PID_Angle_Pair> sorted_intersections;
            if (theta > 0) { // ccw rotation from theta1 to theta2
                for (const auto &pid_angle : edge_intersection_list[ei]) {
                    sorted_intersections.emplace_back(pid_angle.first, compute_rotation_angle(theta1, pid_angle.second));
                }
                std::sort(sorted_intersections.begin(), sorted_intersections.end(), PID_angle_less_than);
            } else { // cw rotation from theta1 to theta2
                for (const auto &pid_angle : edge_intersection_list[ei]) {
                    sorted_intersections.emplace_back(pid_angle.first,
                                                      theta + compute_rotation_angle(theta2, pid_angle.second));
                }
                std::sort(sorted_intersections.begin(), sorted_intersections.end(), PID_angle_greater_than);
            }

            sorted_intersections.emplace_back(edges[ei].id2, theta);
            size_t last_vId = edges[ei].id1;
            double last_angle = 0;
            for (const auto & pid_angle : sorted_intersections) {
                pEdges.emplace_back(SubArc_Edge{last_vId, pid_angle.first, ei,
                                                Circular_Arc(pts[last_vId],pts[pid_angle.first],
                                                             pid_angle.second - last_angle,
                                                             edges[ei].arc.get_center(), edges[ei].arc.get_radius(),
                                                             theta1 + last_angle,
                                                             theta1 + pid_angle.second)});

                last_vId = pid_angle.first;
                last_angle = pid_angle.second;
            }
        }
    }

    //done
}

void Arrangement::decompose_into_cells(const std::vector<Point> &vertices, const std::vector<SubArc_Edge> &edges,
                                       std::vector<std::vector<size_t>> &cells)
{
    // find in-coming and out-going edges for each vertex
    std::vector<std::vector<size_t>> eIn(vertices.size());
    std::vector<std::vector<size_t>> eOut(vertices.size());

    for (int ei = 0; ei < edges.size(); ++ei) {
        eIn[edges[ei].id2].push_back(ei);
        eOut[edges[ei].id1].push_back(ei);
    }

    // find successor for each half-edge
    std::vector<size_t> next_hEdge(edges.size()*2, 0);

    struct Incident_Edge_Data {
        // edge index
        size_t id;
        // is the edge an in-coming edge for the vertex?
        bool is_in;
        // incident angle: angle of the vector tangent to the arc at the vertex
        double angle;
    };
    auto Incident_Edge_Data_greater_than = [](const Incident_Edge_Data &left, const Incident_Edge_Data &right)
    { return left.angle > right.angle; };

    for (int i = 0; i < vertices.size(); ++i) {
        const auto &in_edges = eIn[i];
        const auto &out_edges = eOut[i];
//        if (in_edges.empty()) continue;
        assert(!(in_edges.empty()));

        if (in_edges.size() == 1) {
            // degree-2 vertex
            next_hEdge[2*in_edges[0]+1] = 2*out_edges[0] + 1;
            next_hEdge[2*out_edges[0]] = 2*in_edges[0];
        } else {
            // intersection vertex
            std::vector<Incident_Edge_Data> incident_edge_list;
            for (const auto id : in_edges) {
                double arc_angle = edges[id].arc.get_arc_angle();
                double angle2 = edges[id].arc.get_end_angle();
                double incident_angle = angle2 + ((arc_angle > 0) ? (-M_PI_2) : M_PI_2);
                incident_angle = angle_mod_2PI(incident_angle);
                incident_edge_list.emplace_back(Incident_Edge_Data{id, true, incident_angle});
            }
            for (const auto id : out_edges) {
                double arc_angle = edges[id].arc.get_arc_angle();
                double angle1 = edges[id].arc.get_start_angle();
                double incident_angle = angle1 + ((arc_angle > 0) ? M_PI_2 : (-M_PI_2));
                incident_angle = angle_mod_2PI(incident_angle);
                incident_edge_list.emplace_back(Incident_Edge_Data{id, false, incident_angle});
            }

            // sort incident edges by their incident angle (clockwise order)
            std::sort(incident_edge_list.begin(), incident_edge_list.end(), Incident_Edge_Data_greater_than);

            // find successor for each incident half-edge
            size_t n_incident = incident_edge_list.size();
            for (int j = 0; j < n_incident; ++j) {
                const auto &e = incident_edge_list[j];
                const auto &f = incident_edge_list[(j+1)%n_incident];
                size_t he = e.is_in ? (2*e.id+1) : 2*e.id;
                size_t hf = f.is_in ? 2*f.id : (2*f.id+1);
                next_hEdge[he] = hf;
            }
        }
    }

    // trace half-edges into chains (each chain bounds an arrangement cell)
    trace_chains(next_hEdge, cells);

}

void Arrangement::trace_chains(const std::vector<size_t> &next_hEdge, std::vector<std::vector<size_t>> &chains)
{
    size_t n_hEdge = next_hEdge.size();
    std::vector<bool> is_visited(n_hEdge, false);

    // trace chains
    chains.clear();
    for (size_t i = 0; i < n_hEdge; ++i) {
        if (!is_visited[i]) {
            is_visited[i] = true;
            chains.emplace_back();
            auto &chain = chains.back();
            chain.push_back(i);
            auto next = next_hEdge[i];
            while (!is_visited[next]) {
                is_visited[next] = true;
                chain.push_back(next);
                next = next_hEdge[next];
            }
        }
    }
}