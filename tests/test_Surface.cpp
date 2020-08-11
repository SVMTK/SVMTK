#include <catch.hpp>

#include "Surface.h" 
#include "Slice.h"

TEST_CASE("Points inside and outside surface")
{
    typedef Surface::Point_3 Point_3;
    Surface surface; 
    surface.make_cube(0.,0.,0.,2.0,2.0,2.0,1); 
    Surface other; 
    other.make_cube(-1.,-1.,-1.,1.0,1.0,1.0,1); 
    auto a = surface.vertices_inside(other);
    auto b = surface.vertices_outside(other);
    REQUIRE( a.size()==1);
    REQUIRE( b.size()==7);
    REQUIRE( surface.is_point_inside(Point_3(2.1,2.1,2.1))==false );
}


TEST_CASE("Vertices, points and vectors")
{
    Surface surface; 
    surface.make_cube(0.,0.,0.,2.0,2.0,2.0,2); 
 
    auto vertices = surface.get_vertices();
    REQUIRE( vertices.size()==26);
    auto points = surface.get_points();
    REQUIRE( points.size()==26);

}


TEST_CASE("Region operations")
{
    typedef Surface::Point_3 Point_3;
    typedef Surface::point_vector point_vector;
    point_vector point;
    Surface surface; 
    surface.make_cube(0.,0.,0.,2.0,2.0,2.0,2);  // tetraheder ?? 

    auto vertices = surface.closest_vertices(Point_3(2.0,2.0,2.1), 1);
    REQUIRE( vertices.size()==1);
    point = surface.get_points(vertices);
 
    REQUIRE( (point[0].x()==2.0 and point[0].y()==2.0 and point[0].z()==2.0) );
    surface.adjust_vertices_in_region( vertices.begin() ,vertices.end(), -0.3);
    point = surface.get_points(vertices);
    
    REQUIRE( (point[0].x()==1.9 and point[0].y()==1.8 and point[0].z()==1.8) );

    surface.smooth_laplacian_region( vertices.begin(), vertices.end() , -0.9);
    point = surface.get_points(vertices) ;

    std::cout<<"P3  "<<point[0].x() << " " << point[0].y() << " " << point[0].z()  << std::endl;
    REQUIRE( ( (point[0].x()-2.35)<0.01 and (point[0].y()-1.98)<0.01 and (point[0].z()-1.98)<0.01) );


    surface.smooth_taubin_region( vertices.begin(), vertices.end()      ,30);
    point = surface.get_points(vertices) ;

    REQUIRE( (point[0].x()==1.4 and point[0].y()==1.6 and point[0].z()==1.6) );

}

