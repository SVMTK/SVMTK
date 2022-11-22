#include <catch.hpp>

#include "Surface.h" 

TEST_CASE("Points inside and outside surface")
{
    typedef Surface::Point_3 Point_3;
    Surface surface; 
    surface.make_cube(0.,0.,0.,2.0,2.0,2.0,2.0); 
    Surface other; 
    other.make_cube(-1.,-1.,-1.,1.0,1.0,1.0,2.0); 
    auto a = surface.get_vertices_inside(other);
    auto b = surface.get_vertices_outside(other);

    REQUIRE( a.size()==4);
    REQUIRE( b.size()==10);
    REQUIRE( surface.is_point_inside(Point_3(2.1,2.1,2.1))==false );
}


TEST_CASE("Vertices, points and vectors")
{
    Surface surface; 
    surface.make_cube(0.,0.,0.,2.0,2.0,2.0,1.); 
    auto vertices = surface.get_vertices();
    REQUIRE( vertices.size()==48);
    auto points = surface.get_points();
    REQUIRE( points.size()==48);

}


TEST_CASE("Region operations")
{
    typedef Surface::Point_3 Point_3;
    typedef Surface::point_vector point_vector;
    point_vector point;
    Surface surface; 
    surface.make_cube(0.,0.,0.,2.0,2.0,2.0,2.);  // tetraheder ?? 

    auto vertices = surface.get_closest_vertices(Point_3(2.0,2.0,2.1), 1);
    REQUIRE( vertices.size()==1);
    point = surface.get_points(vertices);
 
    REQUIRE( point[0].x()==Approx(2).margin(1e-12) );
    REQUIRE( point[0].y()==Approx(2).margin(1e-12) );
    REQUIRE( point[0].z()==Approx(2).margin(1e-12) );
 
    surface.adjust_vertices_in_region( vertices.begin() ,vertices.end(), -0.3);
    point = surface.get_points(vertices);
    
    REQUIRE( point[0].x()==Approx(1.826).margin(1e-3) );
    REQUIRE( point[0].y()==Approx(1.826).margin(1e-3) );
    REQUIRE( point[0].z()==Approx(1.826).margin(1e-3) );

    surface.smooth_laplacian_region( vertices.begin(), vertices.end() , -0.9);
    point = surface.get_points(vertices) ;

    REQUIRE( point[0].x()==Approx(2.27).margin(1e-3) );
    REQUIRE( point[0].y()==Approx(2.27).margin(1e-3) );
    REQUIRE( point[0].z()==Approx(2.27).margin(1e-3) );

    surface.smooth_taubin_region( vertices.begin(), vertices.end()      ,50);
    point = surface.get_points(vertices) ;
    REQUIRE( point[0].x()==Approx(1.333).margin(1e-3) );
    REQUIRE( point[0].y()==Approx(1.333).margin(1e-3) );
    REQUIRE( point[0].z()==Approx(1.333).margin(1e-3) );
}




