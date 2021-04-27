//
// Created by Charles Du on 4/26/21.
//

#include <algorithm>
#include "Arrangement.h"

void Arrangement::compute_arrangement(const std::vector<Point> &vertices, const std::vector<Arc_Edge> &edges)
{
    // step 1: subdivide input arcs by their intersections
    std::vector<Point> pts;
    std::vector<SubArc_Edge> pEdges;
    std::vector<bool> is_intersection_point;
    subdivide_polyArc_by_intersection(vertices, edges, pts, pEdges,is_intersection_point);

    // todo: step 2: get all arrangement cells
//    decompose_into_cells(pts, edges, ...)

    // step 3: compute winding numbers for each cell

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
            pEdges.emplace_back(SubArc_Edge{edges[ei].id1, edges[ei].id2, edges[ei].arc.get_arc_angle(), ei});
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
                pEdges.emplace_back(SubArc_Edge{last_vId, pid_angle.first, pid_angle.second - last_angle, ei});
                last_vId = pid_angle.first;
                last_angle = pid_angle.second;
            }
        }
    }

    //done
}