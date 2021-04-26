//
// Created by Charles Du on 4/26/21.
//

#include "Circular_Arc.h"
#include <iostream>

int main(int argc, char const *argv[])
{
    Circular_Arc arc1(Point(0.5, 0.4), Point(-0.5, -0.3), 2.1);
    Circular_Arc arc2(Point(-0.6, 0.5), Point(0.3, 0.5), 2.1);

    std::vector<Intersection_Point> result;
    Circular_Arc::compute_intersection(arc1, arc2, result);

    //
    for (const auto & p : result) {
        std::cout << "------" << std::endl;
        std::cout << "p: (" << p.location.x() << ", " << p.location.y() << ")" << std::endl;
        std::cout << "theta: " << p.angle1 << std::endl;
        std::cout << "phi  : " << p.angle2 << std::endl;
    }

    return 0;
}