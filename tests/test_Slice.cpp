#include <catch.hpp>  // for Approx, operator==, AssertionHandler, SourceLin...
#include <vector>     // for vector
#include "Slice.h"    // for Slice, Slice::Point_2


TEST_CASE("Minimum bounding circle")
{

    typedef Slice::Point_2 Point_2;
    std::vector<Point_2> points; 
    points.push_back(Point_2( 1, 0));  
    points.push_back(Point_2( 0, 1)); 
    points.push_back(Point_2(-1, 0)); 
    points.push_back(Point_2( 0,-1));
    Slice slice;
    slice.add_constraint(points);  
    REQUIRE( slice.get_bounding_circle_radius()==Approx(1.0).margin(1e-3)) ; // less than 4 larger than 2 ?
}

