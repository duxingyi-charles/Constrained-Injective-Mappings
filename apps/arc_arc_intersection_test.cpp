//
// Created by Charles Du on 4/26/21.
//

#include "Circular_Arc.h"
#include <iostream>

int main(int argc, char const *argv[])
{
    Circular_Arc arc1(Point(0.538, 0.32000000000000006), Point(-0.8910065241883679, -0.45399049973954675), -3.2324853071795867);
    Circular_Arc arc2(Point(0.1140000000000001, 0.10999999999999988), Point(0.15643446504023087, -0.9876883405951378), -3.2324853071795867);

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