#ifndef __SURFACE_MESHER_H

#define __SURFACE_MESHER_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> 
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/make_surface_mesh.h>


#include <CGAL/Surface_mesh_traits_generator_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Surface_mesh_triangulation_generator_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>


#include <cstdlib>


#include <CGAL/Surface_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/property_map.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/distance.h>

#include <vector>
#include <fstream>
#include <fstream>
#include <math.h>
#include <functional>


template <typename FT, typename P>
class FT_to_point_function_wrapper : public std::unary_function<P, FT>
{
  //typedef FT (*Implicit_function)(FT, FT, FT);
  typedef std::function<double(double,double,double)> Implicit_function;
  Implicit_function function;
public:
  typedef P Point;
  FT_to_point_function_wrapper(Implicit_function f) : function(f) {}
  FT operator()(Point p) const { return function(p.x(), p.y(), p.z()); }
  ~FT_to_point_function_wrapper(){};
  FT_to_point_function_wrapper();
};


template<typename Mesh , typename Implicit_Function> // Implicit function -> return coordiantes x0,y0,z0-> make center 
void surface_mesher(Mesh& mesh, Implicit_Function func, double& x0 ,double& y0, double& z0 , double bounding_sphere_radius, double angular_bound, double radius_bound, double distance_bound ) 
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Surface_mesh_triangulation_generator_3<Kernel>::Type Tr;
    typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
    typedef Tr::Geom_traits GT;
    typedef Kernel::Point_3 Point_3;
    typedef Kernel::FT FT;

    typedef FT_to_point_function_wrapper<FT, Point_3> Function;
    typedef CGAL::Implicit_surface_3<Kernel, Function> Surface_3;

    // no known conversion for argument 1 from ‘CGAL::Point_3<CGAL::Epeck>’ to ‘const Point_3_& {aka const CGAL::Point_3<CGAL::Epick>&}’
    Tr tr;
    C2t3 c2t3 (tr);
    Function wrapper(func);
    GT::Point_3 bounding_sphere_center(x0, y0, z0);
    GT::FT bounding_sphere_squared_radius = bounding_sphere_radius*bounding_sphere_radius*2;
    GT::Sphere_3 bounding_sphere(bounding_sphere_center,bounding_sphere_squared_radius);

    Surface_3 surface(wrapper,bounding_sphere,1.0e-5 ); 

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angular_bound,radius_bound,distance_bound);

    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3,mesh);

}
template<typename Mesh , typename Implicit_Function>
void surface_mesher(Mesh& mesh, Implicit_Function func,  double bounding_sphere_radius, double angular_bound, double radius_bound, double distance_bound )
{
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef CGAL::Surface_mesh_triangulation_generator_3<Kernel>::Type Tr;

    typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
    typedef Tr::Geom_traits GT;
    typedef Kernel::Point_3 Point_3;
    typedef Kernel::FT FT;
    //Surface_mesh_traits().construct_initial_points_object()(surface_of_sphere_2, CGAL::inserter(tr_3), initial_number_of_points)


    typedef FT_to_point_function_wrapper<FT, Point_3> Function;
    typedef CGAL::Implicit_surface_3<Kernel, Function> Surface_3;

    Tr tr;
    C2t3 c2t3 (tr);
    Function wrapper(func);
    GT::Point_3 bounding_sphere_center(CGAL::ORIGIN);
    GT::FT bounding_sphere_squared_radius = bounding_sphere_radius*bounding_sphere_radius*2;
    GT::Sphere_3 bounding_sphere(bounding_sphere_center,bounding_sphere_squared_radius);
    
    Surface_3 surface(wrapper,bounding_sphere,1.0e-5 ); 

    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angular_bound,radius_bound,distance_bound);


//    CGAL::make_surface_mesh(c2t3, surface, criteria,CGAL::Manifold_with_boundary_tag());
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3,mesh);

}
template<typename Surface>
void poisson_reconstruction(Surface &surface,double angular_bound, double radius_bound, double distance_bound)
{ 
     typedef typename Surface::Kernel Kernel;
     typedef typename Surface::FT FT;
     typedef typename Surface::Point_3 Point_3;
     typedef typename Surface::Vector_3 Vector_3;
     typedef typename Surface::Sphere_3 Sphere_3;
     typedef typename CGAL::Surface_mesh_triangulation_generator_3<Kernel>::Type Tr;
     typedef typename CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
     typedef typename CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function > Implicit_surface_3;
     typedef std::pair<Point_3, Vector_3> Point_with_normal;
     typedef std::vector<Point_with_normal> PointList;
     typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
     typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;


     FT sm_distance(distance_bound) ;
     typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;


     PointList points;
     Tr tr;
     C2t3 c2t3 (tr);
     points = surface.get_points_with_normal();
     Poisson_reconstruction_function function(points.begin(), points.end(), Point_map(), Normal_map());

     if ( ! function.compute_implicit_function() )
     {
        std::cout << "Could not compute implicit function" << std::endl;
        return;
     }
     FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (points, 6 /* knn = 1 ring */,
       CGAL::parameters::point_map (Point_map()));

     Point_3 inner_point = function.get_inner_point();
     Sphere_3 bsphere = function.bounding_sphere();
     FT radius = std::sqrt(bsphere.squared_radius());

     FT sm_sphere_radius = 5.0 * radius;
     FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; 
     Implicit_surface_3 implicit_surface(function,
                      Sphere_3(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error/sm_sphere_radius);
  

     CGAL::Surface_mesh_default_criteria_3<Tr> criteria(angular_bound,radius_bound,distance_bound);


//    CGAL::make_surface_mesh(c2t3, surface, criteria,CGAL::Manifold_with_boundary_tag());
     CGAL::make_surface_mesh(c2t3, implicit_surface, criteria, CGAL::Manifold_tag());
     CGAL::facets_in_complex_2_to_triangle_mesh(c2t3,surface.get_mesh());
 
}
#endif
