//
// Created by Charles Du on 4/21/21.
//

#ifndef TLC_CIRCULAR_ARC_H
#define TLC_CIRCULAR_ARC_H

#include <utility>
#include <vector>

#include "geo_util.h"
#include "Rectangle.h"

template <typename Scalar>
struct Intersection_Point {
    Point<Scalar> location;
    // circular angle in the first arc
    Scalar angle1;
    // circular angle in the second arc
    Scalar angle2;
};

enum Point_Arc_Location { Start, End, Middle };

template <typename Scalar>
class Circular_Arc {
public:
    Circular_Arc(Point<Scalar> p1, Point<Scalar> p2, Scalar theta)
            : start_point(std::move(p1)), end_point(std::move(p2)), arc_angle(theta) { update_arc(); };

    // data constructor
    // user is responsible for input data consistency.
    Circular_Arc(Point<Scalar> p1, Point<Scalar> p2, Scalar theta,
                 Point<Scalar> o, Scalar r, Scalar theta1, Scalar theta2)
                 : start_point(std::move(p1)), end_point(std::move(p2)), arc_angle(theta), center(std::move(o)), radius(r),
                 start_angle(theta1), end_angle(theta2) {};

    ~Circular_Arc() = default;


    Rectangle<Scalar> get_bounding_box() const;

    // compute the most left point on the arc
    // return: pair of <Point p, Point_Arc_Location s>
    std::pair<Point<Scalar>,Point_Arc_Location> get_most_left_point() const;

    // compute the vector tangent to the arc at end point
    Vec2D<Scalar> get_out_tangent_vector() const;
    // compute the vector tangent to the arc at start point
    Vec2D<Scalar> get_in_tangent_vector() const;

    // compute the reverse arc (exchange start and end)
    static Circular_Arc reverse(const Circular_Arc &c);

    // compute the signed area of arc segment
    Scalar get_segment_area() const;

    // compute intersections between arc1 and arc2, save intersection points to result
    static void compute_intersection(const Circular_Arc &arc1, const Circular_Arc &arc2,
                                     std::vector<Intersection_Point<Scalar>>& result);

    // add intersections between arc1 and arc2 to result
    // assume arc1 and arc2 don't share any common end points.
    static void find_all_intersections(const Circular_Arc &arc1, const Circular_Arc &arc2,
                                       std::vector<Intersection_Point<Scalar>>& result);

    // given p is an intersection between arc1 and arc2, add the other intersection (if any) to result
    static void find_other_intersection(const Circular_Arc &arc1, const Circular_Arc &arc2, const Point<Scalar>& p,
                                        std::vector<Intersection_Point<Scalar>>& result);




    Point<Scalar>  get_start_point() const { return start_point; }
    Point<Scalar>  get_end_point() const { return end_point; }
    Scalar get_arc_angle() const { return arc_angle; }
    Scalar get_start_angle() const { return start_angle; }
    Scalar get_end_angle() const { return end_angle; }
    Point<Scalar>  get_center() const { return center; }
    Scalar get_radius() const { return radius; }

private:
    // compute arc center, radius, start and end angle
    void update_arc();


private:
    Point<Scalar> start_point;
    Point<Scalar> end_point;
    Scalar arc_angle;

    Point<Scalar> center;
    Scalar radius{};
    Scalar start_angle{};
    Scalar end_angle{};

};


#endif //TLC_CIRCULAR_ARC_H
