#ifndef CGALSurface_H


#define CGALSurface_H

#ifndef BOOST_PARAMETER_MAX_ARITY
# define BOOST_PARAMETER_MAX_ARITY 12
#endif

class Implicit_function;

#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <memory>
#include <array>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Mesh_polyhedron_3.h>






class CGALSurface
{
   public:



    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel; // Change to exact -> copy face_ graph to inexact ??
    typedef Kernel::Point_3 Point_3;
    typedef CGAL::Surface_mesh<Point_3> Mesh;
    typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
    typedef Tr::Geom_traits GT;
    typedef Kernel::Vector_3 Vector_3;




    typedef boost::graph_traits<Mesh>::halfedge_descriptor    halfedge_descriptor;
    typedef boost::graph_traits<Mesh>::face_descriptor        face_descriptor;
    typedef boost::graph_traits<Mesh>::vertex_descriptor      vertex_descriptor;
    typedef std::vector<vertex_descriptor>                    vertex_vector;
    typedef CGAL::Side_of_triangle_mesh<Mesh,Kernel> Inside; //


    CGALSurface(const std::string filename);

    CGALSurface( Implicit_function implicit_function,
                 double bounding_sphere_radius,
                 double angular_bound,
                 double radius_bound,
                 double distance_bound);


        void operator^=( CGALSurface& other);

        void operator+=( CGALSurface& other );

        void operator-=( CGALSurface& other );

        int fill_holes();

        bool triangulate_faces();
      
        void stitch_borders();
         
        void insert_surface(CGALSurface& surface);

        void isotropic_remeshing(double target_edge_length, unsigned int nb_iter, bool protect_border);
       
        void adjust_boundary(const double c);

        void smooth_laplacian(const double c);

        template<typename InputIterator>
        void adjusting_boundary_region(InputIterator begin , InputIterator  end, const double c);

        template<typename InputIterator>
        void smooth_laplacian_region(InputIterator begin , InputIterator end ,const int c);

        vertex_vector points_inside(CGALSurface& other);

        vertex_vector points_outside(CGALSurface& other);

        Mesh& get_mesh();
        //Polyhedron& polyhedron();

        bool self_intersections();

        template< typename Polyhedron_3>
        void get_polyhedron(Polyhedron_3 &polyhedron_3);

        void save(const char* outpath);

        void preprocess(double target_edge_length, int nb_iter);


        void fair(CGALSurface::vertex_vector vector  );

        template<typename Implicit_function>
        CGALSurface(Implicit_function implicit_function,
             double bounding_sphere_radius,
             double angular_bound=30.,
             double radius_bound=0.1 ,
             double distance_bound=0.1   );

        ~CGALSurface(){}
   private :
       Mesh mesh;


};

#endif
