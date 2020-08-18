#include <catch.hpp>

#include "Slice.h" 



TEST_CASE("Minimum bounding circle")
{

    typedef Slice::Point_2 Point_2;
    std::vector<Point_2> points; 
    points.push_back(Point_2( 1, 0));  
    points.push_back(Point_2( 0, 1)); 
    points.push_back(Point_2(-1, 0)); 
    points.push_back(Point_2( 0,-1));
    Slice slice;
    slice.add_constraints(points);  
    REQUIRE( slice.get_bounding_circle_radius()==Approx(1).margin(1e-3)) ; // less than 4 larger than 2 ?
}

