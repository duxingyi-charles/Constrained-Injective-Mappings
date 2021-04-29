//
// Created by Charles Du on 4/21/21.
//

#include "Circular_Arc.h"

void Circular_Arc::update_arc() {
    VarVec2d vec = end_point - start_point;
    radius = vec.norm()/(2*sin(abs(arc_angle)/2));

    VarVec2d perpVec(-vec.y(), vec.x());
    center = (start_point + end_point)/2 + perpVec/(2*tan(arc_angle/2));

    start_angle = compute_vector_angle(start_point - center);
    end_angle = start_angle + arc_angle;
}


Rectangle Circular_Arc::get_bounding_box() const {
    var xmin, xmax, ymin, ymax;
    // left
    if (is_angle_between(M_PI, start_angle, end_angle)) {
        xmin = center.x() - radius;
    } else {
//        xmin = fmin(start_point.x(), end_point.x());
        xmin = start_point.x() < end_point.x() ? start_point.x() : end_point.x();
    }
    // right
    if (is_angle_between(0, start_angle, end_angle)) {
        xmax = center.x() + radius;
    } else {
//        xmax = fmax(start_point.x(), end_point.x());
        xmax = start_point.x() > end_point.x() ? start_point.x() : end_point.x();
    }
    // bottom
    if (is_angle_between(3*M_PI/2, start_angle, end_angle)) {
        ymin = center.y() - radius;
    } else {
//        ymin = fmin(start_point.y(), end_point.y());
        ymin = start_point.y() < end_point.y() ? start_point.y() : end_point.y();
    }
    // top
    if (is_angle_between(M_PI/2, start_angle, end_angle)) {
        ymax = center.y() + radius;
    } else {
//        ymax = fmax(start_point.y(), end_point.y());
        ymax = start_point.y() > end_point.y() ? start_point.x() : end_point.y();
    }
    //
    return Rectangle(Point(xmin,ymin), Point(xmax,ymax));
}

std::pair<Point,Point_Arc_Location> Circular_Arc::get_most_left_point() const
{
    if (is_angle_between(M_PI, start_angle, end_angle)) {
        return std::make_pair(Point(center.x()-radius,center.y()), Middle);
    } else {
        if (start_point.x() < end_point.x()) {
            return std::make_pair(start_point, Start);
        } else {
            return std::make_pair(end_point, End);
        }
    }
}

VarVec2d Circular_Arc::get_out_tangent_vector() const
{
    if (arc_angle > 0) {
        return rotate_90deg(end_point - center);
    } else {
        return rotate_90deg(center - end_point);
    }
}

VarVec2d Circular_Arc::get_in_tangent_vector() const
{
    if (arc_angle > 0) {
        return rotate_90deg(start_point - center);
    } else {
        return rotate_90deg(center - start_point);
    }
}

Circular_Arc Circular_Arc::reverse(const Circular_Arc &c) {
    return Circular_Arc(c.end_point, c.start_point, -c.arc_angle,
                        c.center, c.radius, c.end_angle, c.start_angle);
}

void Circular_Arc::compute_intersection(const Circular_Arc &arc1, const Circular_Arc &arc2,
                          std::vector<Intersection_Point>& result)
{
    result.clear();
    const Point &p1 = arc1.start_point;
    const Point &p2 = arc1.end_point;
    const Point &q1 = arc2.start_point;
    const Point &q2 = arc2.end_point;

    if (p1==q1 && p2==q2) {
        if (arc1.arc_angle == arc2.arc_angle) {
            // two arcs completely overlap, no need to create intersection
            return;
        } else {
            result.emplace_back(Intersection_Point{p1, arc1.start_angle, arc2.start_angle});
        }
    } else if (p1==q2 && q1==p2) {
        return;
    } else if (p1 == q1) {
        result.emplace_back(Intersection_Point{p1, arc1.start_angle, arc2.start_angle});
        find_other_intersection(arc1, arc2, p1, result);
    } else if (p1 == q2) {
        find_other_intersection(arc1, arc2, p1, result);
    } else if (p2 == q2 || p2 == q1) {
        find_other_intersection(arc1, arc2, p2, result);
    } else { // arc1 and arc2 don't share common end points
        find_all_intersections(arc1, arc2, result);
    }
}

void Circular_Arc::find_all_intersections(const Circular_Arc &arc1, const Circular_Arc &arc2,
                            std::vector<Intersection_Point>& result)
{
    const Point &o1 = arc1.center;
    const Point &o2 = arc2.center;
    var r1 = arc1.radius;
    var r2 = arc2.radius;

    var dist_o1o2 = (o1-o2).norm();
    if (dist_o1o2 > r1 + r2 || dist_o1o2 < abs(r1-r2) || dist_o1o2 == 0) {
        // when dist_o1o2==0, the two arcs could have some overlap, but we don't treat that as intersection
        return;
    }

    // compute intersection of two circles
    var c = 2 * r1 * dist_o1o2 * dist_o1o2;
    var o1o2_diff_x = o1.x() - o2.x();
    var o1o2_diff_y = o1.y() - o2.y();
    var diff2 = r2*r2 - r1*r1 - dist_o1o2*dist_o1o2;
    var sqroot = sqrt((r1+dist_o1o2+r2)*(r1+dist_o1o2-r2)*(r2+r1-dist_o1o2)*(r2-r1+dist_o1o2));

    std::vector<var> sinTheta_list, cosTheta_list; // theta indicates circular angles in arc1
    if (abs(o1.x()-o2.x()) > abs(o1.y()-o2.y())) {
        var a = diff2 * o1o2_diff_y;
        var b = o1o2_diff_x * sqroot;
        var sin1 = (a+b)/c;
        var sin2 = (a-b)/c;
        if (abs(sin1) <= 1) {
            sinTheta_list.push_back(sin1);
        }
        if (abs(sin2) <= 1) {
            sinTheta_list.push_back(sin2);
        }

        for (const var &s : sinTheta_list) {
            cosTheta_list.emplace_back((diff2 - 2*r1*o1o2_diff_y * s)/(2*r1*o1o2_diff_x));
        }
    } else {
        var a = diff2 * o1o2_diff_x;
        var b = o1o2_diff_y * sqroot;
        var cos1 = (a+b)/c;
        var cos2 = (a-b)/c;
        if (abs(cos1) <= 1) {
            cosTheta_list.push_back(cos1);
        }
        if (abs(cos2) <= 1) {
            cosTheta_list.push_back(cos2);
        }

        for (const var &s : cosTheta_list) {
            sinTheta_list.emplace_back((diff2 - 2*r1*o1o2_diff_x * s)/(2*r1*o1o2_diff_y));
        }
    }

    std::vector<var> sinPhi_list, cosPhi_list;  //phi indicates circular angles in arc2
    for (const var &s : sinTheta_list) {
        sinPhi_list.emplace_back((o1o2_diff_y + r1 * s)/r2);
    }
    for (const var &s : cosTheta_list) {
        cosPhi_list.emplace_back((o1o2_diff_x + r1 * s)/r2);
    }

    // extract circular angles from sin and cos
    std::vector<var> theta_list, phi_list;
    for (int i = 0; i < cosTheta_list.size(); ++i) {
        theta_list.push_back(compute_vector_angle(Point(cosTheta_list[i], sinTheta_list[i])));
    }
    for (int i = 0; i < cosPhi_list.size(); ++i) {
        phi_list.push_back(compute_vector_angle(Point(cosPhi_list[i], sinPhi_list[i])));
    }

    // select circular angles in the range of input arc angles
    std::vector<var> final_theta_list, final_phi_list;
    for (int i = 0; i < theta_list.size(); ++i) {
        if (is_angle_between(theta_list[i], arc1.start_angle, arc1.end_angle) &&
            is_angle_between(phi_list[i], arc2.start_angle, arc2.end_angle)) {
            final_theta_list.push_back(theta_list[i]);
            final_phi_list.push_back(phi_list[i]);
        }
    }
    
    // record intersections
    for (int i = 0; i < final_theta_list.size(); ++i) {
        var theta = final_theta_list[i];
        var phi = final_phi_list[i];
        VarVec2d vec1(cos(theta), sin(theta));
        VarVec2d vec2(cos(phi), sin(phi));
        result.emplace_back(Intersection_Point{
                (o1 + r1 * vec1 + o2 + r2 * vec2)/2,
                final_theta_list[i],
                final_phi_list[i]});
    }

}

void Circular_Arc::find_other_intersection(const Circular_Arc &arc1, const Circular_Arc &arc2, const Point& p,
                             std::vector<Intersection_Point>& result)
{
    const Point &o1 = arc1.center;
    const Point &o2 = arc2.center;

    if (o1 == o2) { // we don't treat overlap as intersection
        return;
    }

    if ((p.x()-o1.x())*(o2.y()-o1.y())-(p.y()-o1.y())*(o2.x()-o1.x()) == 0) {
        // p lies on line o1o2, there is no other intersection
        return;
    }

    // compute the other intersection between the two circles
    var t = ((p-o1).dot(o2-o1))/((o2-o1).dot(o2-o1));
    Point proj_p = o1 + t * (o2 - o1);
    Point q = 2 * proj_p - p;
    var theta = compute_vector_angle(q - o1);
    var phi   = compute_vector_angle(q - o2);

    // if the intersection lies on the arcs, save it to result
    if (is_angle_between(theta, arc1.start_angle, arc1.end_angle)
    &&  is_angle_between(phi  , arc2.start_angle, arc2.end_angle)) {
        result.emplace_back(Intersection_Point{q, theta, phi});
    }
}

var Circular_Arc::get_segment_area() const {
    return 0.5 * radius * radius * (arc_angle - sin(arc_angle));
}
