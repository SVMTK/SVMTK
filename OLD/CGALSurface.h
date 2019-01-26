#ifndef CGALSurface_H
#define CGALSurface_H

#define BOOST_PARAMETER_MAX_ARITY 12

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "utils.h"

#include <iostream>
#include <vector>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

// Isotropic remeshing
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

// Edge collapse -- Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

// copy_face_graph
#include <CGAL/boost/graph/copy_face_graph.h>

// points inside/outside
#include <CGAL/Side_of_triangle_mesh.h>

// Corefine and compute difference
#include <CGAL/Polygon_mesh_processing/corefinement.h>

// Poisson reconstruction
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/poisson_surface_reconstruction.h>


// Poisson reconstruction definitions- NB! Another Kernel



class CGALSurface {
    // What TODO: were these two for?
    /*     typedef CGAL::Surface_mesh_default_triangulation_3 Tr; */
    /*     typedef Tr::Geom_traits GT; */
    public:
        typedef CGAL::Exact_predicates_inexact_constructions_kernel ReconstructKernel;
        typedef CGAL::Polyhedron_3<ReconstructKernel> RPolyhedron;
        typedef ReconstructKernel::Point_3 RPoint;
        typedef ReconstructKernel::Vector_3 RVector;
        typedef std::pair<RPoint, RVector> RPwn;

        typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel; // Change to exact -> copy face_ graph to inexact ??
        typedef Kernel::Point_3 Point_3;
        typedef CGAL::Surface_mesh<Point_3> Mesh;
        typedef Kernel::Vector_3 Vector_3;
        typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
        typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
        typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
        typedef CGAL::Side_of_triangle_mesh<Mesh, Kernel> inside; //  TODO: Cnange name, unintuitive



        CGALSurface(); // empty constructor

        CGALSurface(const std::string f);

        CGALSurface(Polyhedron &);

        /* CGALSurface( Implicit_function implicit_function, */
        /*              double bounding_sphere_radius, */
        /*              double angular_bound, */
        /*              double radius_bound, */
        /*              double distance_bound); */

        /* template<typename Implicit_function> */
        /* CGALSurface(Implicit_function implicit_function, */
        /*          double bounding_sphere_radius, */
        /*          double angular_bound=30., */
        /*          double radius_bound=0.1 , */
        /*          double distance_bound=0.1   ); */

        void load(const std::string f);


        template<typename Polyhedron_3>
        void get_polyhedron(Polyhedron_3 &polyhedron_3);

        void operator^=(CGALSurface &other); 
        void operator+=(CGALSurface &other); 
        void operator-=(CGALSurface &other); 


        void surface_intersection(CGALSurface &other);
        void surface_union(CGALSurface &other);
        void surface_difference(CGALSurface &other);


        Mesh& get_mesh();

        int fill_holes();

        bool triangulate_faces();

        void stitch_borders(); 

        void insert_surface(CGALSurface& surface);    // TODO

        void isotropic_remeshing(const double, const unsigned int, const bool);

        void adjust_boundary(const double);

        void smooth_laplacian(const double);

        void smooth_taubin(const size_t);

        template<typename InputIterator>
        void adjusting_boundary_region(InputIterator, InputIterator, const double);

        template<typename InputIterator>
        void smooth_laplacian_region(InputIterator, InputIterator, const int);

        std::vector<vertex_descriptor> points_inside(CGALSurface &);

        std::vector<vertex_descriptor> points_outside(CGALSurface& other);

        bool self_intersections();

        int num_self_intersections();

        void save(const std::string);

        int collapse_edges(const double stop_ratio);

        void preprocess(const double, const int);

        void fair(std::vector<vertex_descriptor> &);

        // Experimental reconstruction
        void reconstruct_surface(
            const double sm_angle = 20.0,
            const double sm_radius = 30.0,
            const double sm_distance = 0.375);

        int num_faces() const;

        int num_edges() const;

        int num_vertices() const;

        // TODO: What does this do?
        void clear(){ mesh.clear();}

        // TODO: Add this for tomorrow
        void fix_close_junctures(double c); 

        // TODO
        void split_edges(double  target_edge_length); 
        void make_cylinder( double x0, double y0, double  z0,  double x1, double y1, double z1, double radius,  int number_of_segments=360) ; 
        void make_cone( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0 , double r1,  int number_of_segments=360) ; 
        void make_cube( double x0, double y0, double  z0,  double x1, double y1, double z1); 
        void make_sphere( double x0, double y0, double  z0, double r0); 
        void insert_mesh(CGALSurface& surf){mesh+=surf.get_mesh();} 
        void insert_points(std::vector<Point_3>& points) ; 

        template<typename Implicit_function> 
        CGALSurface(Implicit_function implicit_function, 
              double bounding_sphere_radius, 
              double angular_bound=30., 
              double radius_bound=0.1 , 
              double distance_bound=0.1   ); 

        ~CGALSurface(){}
    private:
        Mesh mesh;
};



#endif
