#ifndef __SURFACE_MESHER_H

#define __SURFACE_MESHER_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> 
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_traits_generator_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Surface_mesh_triangulation_generator_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
//#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Timer.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/compute_average_spacing.h>

#include <deque>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <functional>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// needed

#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/boost/graph/copy_face_graph.h>


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
    typedef Kernel::Sphere_3 Sphere_3;
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
    typedef Kernel::Sphere_3 Sphere_3;
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

struct Counter {
    std::size_t i, N;
    Counter(std::size_t N) : i(0), N(N) {}

    void operator()() {
        i++;
    }
};


struct InsertVisitor {
    Counter& c;
    InsertVisitor(Counter& c) : c(c) {}

    void before_insertion() {
        c();
    }
};


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

template < typename Mesh>
void poisson_reconstruction(Mesh &mesh,
        const double sm_angle = 20.0,
        const double sm_radius = 100.0,
        const double sm_distance = 0.25,
        const double approximation_ratio = 0.02,
        const double average_spacing_ratio = 5.0) {



        typedef CGAL::Exact_predicates_inexact_constructions_kernel ReconstructKernel;

        // Simple geometric types
        typedef ReconstructKernel::FT FT;
        typedef ReconstructKernel::Point_3 Point;
        typedef ReconstructKernel::Vector_3 Vector;
        typedef CGAL::Point_with_normal_3<ReconstructKernel> Point_with_normal;
        typedef ReconstructKernel::Sphere_3 Sphere;
        typedef std::deque<Point_with_normal> PointList;

        // polyhedron -- Moved to Surface
        typedef CGAL::Polyhedron_3<ReconstructKernel> RPolyhedron;

        // Poisson implicit function
        typedef CGAL::Poisson_reconstruction_function<ReconstructKernel> Poisson_reconstruction_function;

        // Surface mesher
        typedef CGAL::Surface_mesh_default_triangulation_3 STr;
        typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
        /* typedef CGAL::Poisson_implicit_surface_3<ReconstructKernel, Poisson_reconstruction_function> Surface_3; */
        typedef CGAL::Implicit_surface_3<ReconstructKernel, Poisson_reconstruction_function> Surface_3;







    //***************************************
    // decode parameters
    //***************************************

    // usage
    // Poisson options
    /* FT sm_angle = 20.0; // Min triangle angle (degrees). */
    /* FT sm_radius = 100; // Max triangle size w.r.t. point set average spacing. */
    /* FT sm_distance = 0.25; // Approximation error w.r.t. point set average spacing. */
    /* double sm_angle = 20.0; // Min triangle angle (degrees). */
    /* double sm_radius = 100; // Max triangle size w.r.t. point set average spacing. */
    /* double sm_distance = 0.25; // Approximation error w.r.t. point set average spacing. */
    /* double approximation_ratio = 0.02; */
    /* double average_spacing_ratio = 5; */

    CGAL::Timer task_timer; task_timer.start();
    PointList points;
    for (auto vd : vertices(mesh)) {
        const Point& p = mesh.point(vd);
        Vector n = CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, mesh );
        points.push_back(Point_with_normal(p, n));
    }

    //***************************************
    // Checks requirements
    //***************************************

    std::size_t nb_points = points.size();
    if (nb_points == 0) {
        std::cerr << "Error: empty point set" << std::endl;
    }

    bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if (!points_have_normals) {
        std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
    }

    CGAL::Timer reconstruction_timer; reconstruction_timer.start();

    Counter counter(std::distance(points.begin(), points.end()));
    InsertVisitor visitor(counter) ;

    //***************************************
    // Computes implicit function  error
    //***************************************

    std::cout << "Computes Poisson implicit function...\n";

    // Creates implicit function from the read points.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    Poisson_reconstruction_function function(
            points.begin(), points.end(),
            CGAL::make_identity_property_map(PointList::value_type()),
            CGAL::make_normal_of_point_with_normal_map(PointList::value_type()),
            visitor);

    CGAL::Eigen_solver_traits<Eigen::ConjugateGradient<CGAL::Eigen_sparse_symmetric_matrix<double>::EigenType> > solver;
    bool implicit_success = function.compute_implicit_function(
            solver, visitor,
            approximation_ratio,
            average_spacing_ratio);
    if (!implicit_success) {
        std::cerr << "Error: cannot compute implicit function" << std::endl;
    }

    // Prints status
    std::cerr << "Total implicit function (triangulation+refinement+solver): " << task_timer.time() << " seconds\n";
    task_timer.reset();

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cout << "Surface meshing...\n";

    // Computes average spacing
    FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(
            points, 6 /* knn = 1 ring */);

    // Gets one point inside the implicit surface
    Point inner_point = function.get_inner_point();
    FT inner_point_value = function(inner_point);
    if (inner_point_value >= 0.0) {
        std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
    }

    // Gets implicit function's radius
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());


    // CALL SURFACE MESHER
    // Defines the implicit surface: requires defining a
    // conservative bounding sphere centered at inner point.
    FT sm_sphere_radius = 5.0*radius;
    FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
            Sphere(inner_point, sm_sphere_radius*sm_sphere_radius),
            sm_dichotomy_error/sm_sphere_radius);

    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(
            sm_angle,  // Min triangle angle (degrees)
            sm_radius*average_spacing,  // Max triangle size
            sm_distance*average_spacing); // Approximation error

    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(
            c2t3,                                 // reconstructed mesh
            surface,                              // implicit surface
            criteria,                             // meshing criteria
            CGAL::Manifold_with_boundary_tag());  // require manifold mesh

    // Prints status
    std::cout << "Surface meshing: " << task_timer.time() << " seconds, "
              << tr.number_of_vertices() << " output vertices" << std::endl;
    task_timer.reset();

    if (tr.number_of_vertices() == 0) {
    }


    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3,mesh);

    // Prints total reconstruction duration
    std::cout << "Total reconstruction (implicit function + meshing): " << reconstruction_timer.time() << " seconds\n";
    //return new Surface(tmp_mesh);
}



#endif
