#ifndef CGALSurface_H


#define CGALSurface_H

#ifndef BOOST_PARAMETER_MAX_ARITY
# define BOOST_PARAMETER_MAX_ARITY 12
#endif
class Implicit_function; // forward declartion 
                        

// surface mesher 
// #include <CGAL/make_surface_mesh.h>
// #include <CGAL/Implicit_surface_3.h>
// #include <CGAL/Complex_2_in_triangulation_3.h>
//#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
//#include <array> 



//#include <iterator>
//#include <memory>


//#include <CGAL/Surface_mesh_default_criteria_3.h>


// needed headers
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <vector>
#include <CGAL/Side_of_triangle_mesh.h> // 1 typedef 
#include <CGAL/boost/graph/copy_face_graph.h>
//////////?? 









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
    typedef CGAL::Side_of_triangle_mesh<Mesh,Kernel> Inside; // declear inside implementation ? 
    typedef Mesh::Vertex_index Index;
    typedef std::vector<Point_3>  Polyline_3; 
    typedef std::vector<Polyline_3> Polylines; 
    typedef Kernel::FT FT;

    CGALSurface(); // empty constructor


    CGALSurface(const std::string  filename);

    CGALSurface( Implicit_function implicit_function,
                 double bounding_sphere_radius,
                 double angular_bound,
                 double radius_bound,
                 double distance_bound);

    template<typename Implicit_function>
    CGALSurface(Implicit_function implicit_function,
             double bounding_sphere_radius,
             double angular_bound=30.,
             double radius_bound=0.1 ,
             double distance_bound=0.1   );

    ~CGALSurface(){}


        void split_edges(double  target_edge_length);

        void make_cylinder( double x0, double y0, double  z0,  double x1, double y1, double z1, double radius,  int number_of_segments=360) ;

        void make_cone( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0 , double r1,  int number_of_segments=360) ;

        void make_cube( double x0, double y0, double  z0,  double x1, double y1, double z1);

 
        void make_sphere( double x0, double y0, double  z0, double r0);

        void operator^=( CGALSurface& other);

        void operator+=( CGALSurface& other ); 

        void operator-=( CGALSurface& other );

        void surface_intersection( CGALSurface& other);

        void surface_difference( CGALSurface& other);

        void surface_union( CGALSurface& other);

        void smooth_taubin(const size_t nb_iter); 

        int num_self_intersections();

        int collapse_edges(const double stop_ratio);

        int fill_holes();

        bool triangulate_faces();
      
        void stitch_borders();
         
        void insert_surface(CGALSurface& surface);

        void isotropic_remeshing(double target_edge_length, unsigned int nb_iter, bool protect_border);
       
        void adjust_boundary(const double c);

        void smooth_laplacian(const double c);

        template<typename InputIterator>  // -> operation on surface_mesh-> index-based 
        void adjusting_boundary_region(InputIterator begin , InputIterator  end, const double c);

        template<typename InputIterator>
        void smooth_laplacian_region(InputIterator begin , InputIterator end ,const int c);

        vertex_vector points_inside(CGALSurface& other);

        vertex_vector points_outside(CGALSurface& other);

        Mesh& get_mesh();


        int num_faces() const {return mesh.number_of_faces();}


        int num_edges() const { return mesh.number_of_edges();}


        int num_vertices() const {return mesh.number_of_vertices();}

        bool insert_mesh(CGALSurface& surf){mesh+=surf.get_mesh();}
        
        bool self_intersections();

        template< typename Polyhedron_3>  // Template to address the different kernels in CGAL
        void get_polyhedron(Polyhedron_3 &polyhedron_3 ){CGAL::copy_face_graph(mesh,polyhedron_3);}

        void save(const char* outpath);

        void preprocess(double target_edge_length, int nb_iter);

        void fair(CGALSurface::vertex_vector vector  );

        void poisson_reconstruction() ;

        void clear(){ mesh.clear();}

        void insert_points(std::vector<Point_3>& points) ;

        void fix_close_junctures(double c);

   //private :
   protected:
       Mesh mesh;


};

#endif
