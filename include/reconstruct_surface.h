#ifndef ReconstructSurface_H
#define ReconstructSurface_H
// poisson_reconstruction.cpp

//----------------------------------------------------------
// Poisson Delaunay Reconstruction method.
// Reads a point set or a mesh's set of vertices, reconstructs a surface using Poisson,
// and saves the surface.
// Output format is .off.
//----------------------------------------------------------
// poisson_reconstruction file_in file_out [options]

#include "CGALSurface.h"

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Poisson_implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/compute_average_spacing.h>

/* #include <deque> */
/* #include <cstdlib> */
/* #include <fstream> */
/* #include <math.h> */

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel -- Moved to CGALSurface
/* typedef CGAL::Exact_predicates_inexact_constructions_kernel ReconstructKernel; */

// Simple geometric types
typedef ReconstructKernel::FT FT;
typedef ReconstructKernel::Point_3 Point;
typedef ReconstructKernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<ReconstructKernel> Point_with_normal;
typedef ReconstructKernel::Sphere_3 Sphere;
typedef std::deque<Point_with_normal> PointList;

// polyhedron -- Moved to CGALSurface
typedef CGAL::Polyhedron_3<ReconstructKernel> Polyhedron;

// Poisson implicit function
typedef CGAL::Poisson_reconstruction_function<ReconstructKernel> Poisson_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Poisson_implicit_surface_3<ReconstructKernel, Poisson_reconstruction_function> Surface_3;

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


CGALSurface* reconstruct_surface(Polyhedron &input_mesh) {
    //***************************************
    // decode parameters
    //***************************************

    // usage
    // Poisson options
    FT sm_angle = 20.0; // Min triangle angle (degrees).
    FT sm_radius = 100; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.25; // Approximation error w.r.t. point set average spacing.
    double approximation_ratio = 0.02;
    double average_spacing_ratio = 5;

    CGAL::Timer task_timer; task_timer.start();
    PointList points;
    for (boost::graph_traits<Polyhedron>::vertex_descriptor vd: vertices(input_mesh)) {
        const Point& p = vd->point();
        Vector n = CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, input_mesh);
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
    // Computes implicit function
    //***************************************

    std::cout << "Computes Poisson implicit function...\n";

    // Creates implicit function from the read points.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    Poisson_reconstruction_function function(
            points.begin(), points.end(),
            CGAL::make_identity_property_map(PointList::value_type()),
            CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()),
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
            points.begin(), points.end(), 6 /* knn = 1 ring */);

    // Gets one point inside the implicit surface
    Point inner_point = function.get_inner_point();
    FT inner_point_value = function(inner_point);
    if (inner_point_value >= 0.0) {
        std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
    }

    // Gets implicit function's radius
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());

    // Defines the implicit surface: requires defining a
    // conservative bounding sphere centered at inner point.
    FT sm_sphere_radius = 5.0*radius;
    FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(
            function,
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

    // Converts to polyhedron
    Polyhedron output_mesh;
    CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);

    // Prints total reconstruction duration
    std::cout << "Total reconstruction (implicit function + meshing): " << reconstruction_timer.time() << " seconds\n";
    return new CGALSurface(output_mesh);
}


CGALSurface reconstruct_surface(CGALSurface &surface) {
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);
    return *reconstruct_surface(polyhedron);
}

#endif
