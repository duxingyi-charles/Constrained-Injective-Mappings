//
// Created by Charles Du on 4/21/21.
//

#ifndef TLC_CIRCULAR_ARC_H
#define TLC_CIRCULAR_ARC_H

#include <vector>

#include "geo_util.h"
#include "Rectangle.h"


struct Intersection_Point {
    Point location;
    // circular angle in the first arc
    double angle1;
    // circular angle in the second arc
    double angle2;
};

class Circular_Arc {
public:
    Circular_Arc(const Point &p1, const Point &p2, double theta)
            : start_point(p1), end_point(p2), arc_angle(theta) { update_arc(); };

    ~Circular_Arc() = default;


    Rectangle get_bounding_box();

    // compute intersections between arc1 and arc2, save intersection points to result
    static void compute_intersection(const Circular_Arc &arc1, const Circular_Arc &arc2,
                                     std::vector<Intersection_Point>& result);

    // add intersections between arc1 and arc2 to result
    // assume arc1 and arc2 don't share any common end points.
    static void find_all_intersections(const Circular_Arc &arc1, const Circular_Arc &arc2,
                                       std::vector<Intersection_Point>& result);

    // given p is an intersection between arc1 and arc2, add the other intersection (if any) to result
    static void find_other_intersection(const Circular_Arc &arc1, const Circular_Arc &arc2, const Point& p,
                                        std::vector<Intersection_Point>& result);

private:
    // compute arc center, radius, start and end angle
    void update_arc();


private:
    Point start_point;
    Point end_point;
    double arc_angle;

    Point center;
    double radius;
    double start_angle;
    double end_angle;

};


#endif //TLC_CIRCULAR_ARC_H
