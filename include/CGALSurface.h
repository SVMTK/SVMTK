#ifndef CGALSurface_H
#define CGALSurface_H

#ifndef BOOST_PARAMETER_MAX_ARITY
# define BOOST_PARAMETER_MAX_ARITY 12
#endif

// Local
#include "CGALSlice.h" 
#include "surface_mesher.h"
#include "reconstruct_surface.h"
#include "utils.h"

// BOOST
#include <CGAL/boost/graph/copy_face_graph.h>

//  STL
#include <assert.h>
#include <iterator>
#include <vector>

// TODO: Clean up PMP  header files
//---------------------------------------------------------
//---------------------CLEAN UP-----------------------------------
//---------------------------------------------------------
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/IO/Triangulation_off_ostream_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Optimal_transportation_reconstruction_2.h>

//---------------------------------------------------------
// CGAL- Polygon_mesh_processing
//---------------------------------------------------------
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_slicer.h>

// CGAL accelerated queries
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/squared_distance_3.h> 

#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Poisson_reconstruction_function.h>

// CGAL Surface mesh simplificiation 
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

// ?? 
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/poisson_surface_reconstruction.h>

// CGAL IO
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

// CGAL
/* #include <CGAL/Exact_predicates_exact_constructions_kernel.h> */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>

#include <CGAL/Side_of_triangle_mesh.h>


struct Sphere_wrapper {
    // TODO : Sphere uses a implicit function, and can be extended to other implict fucntion with smooth boundaries.  
    public:
        static double radius ;
        static double x0;
        static double y0;
        static double z0;

        static double function(const double x, const double y, const double z) {
            return (x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0) - radius*radius;
        }
};

double Sphere_wrapper::radius = 0;
double Sphere_wrapper::x0 = 0;
double Sphere_wrapper::y0 = 0;
double Sphere_wrapper::z0 = 0;


template< typename CGALSurface>  // FIXME:
void surface_overlap(CGALSurface& surf1, CGALSurface& surf2)
{
    typedef typename CGALSurface::vertex_vector vertex_vector;
    typedef typename CGALSurface::vertex_descriptor vertex_descriptor;
    // distance from other surface

    std::map< vertex_descriptor, double > map1;
    std::map< vertex_descriptor, double > map2;

    vertex_vector p1 = surf1.get_vertices();
    vertex_vector p2 = surf2.get_vertices();

    p1 = surf2.inside(surf1, p1);
    p2 = surf1.inside(surf2, p2);

    int iter =0;
    while (!p2.empty() or !p1.empty())
    {
         map1 = surf1.shortest_edge_map(p1,-0.1);
         map2 = surf2.shortest_edge_map(p2,-0.1);

         surf1.adjusting_boundary_region(map1.begin() ,map1.end());
         surf2.adjusting_boundary_region(map2.begin() ,map2.end());

         // After adjusting the boundary smoothing is needed, taubin is least volatile.
         surf1.smooth_taubin_region(p1.begin(), p1.end(), 2);
         surf2.smooth_taubin_region(p2.begin(), p2.end(), 2);

         p1 = surf2.inside(surf1, p1);
         p2 = surf1.inside(surf2, p2);

         if (iter++ > 300)
             break;
    }
}


template< typename CGALSurface> 
void surface_overlap(CGALSurface& surf1 , CGALSurface& surf2, CGALSurface& domain)
{
    typedef typename CGALSurface::vertex_vector vertex_vector;
    typedef typename CGALSurface::vertex_descriptor vertex_descriptor;

    std::map< vertex_descriptor, double > map1;
    std::map< vertex_descriptor, double > map2;

    vertex_vector p1 = surf1.get_vertices();
    vertex_vector p2 = surf2.get_vertices();

    p1 = surf2.inside(surf1, p1);
    p2 = surf1.inside(surf2, p2);
    p1 = domain.outside(surf1, p1);
    p2 = domain.outside(surf2, p2);
    int iter = 0;

    while (!p2.empty() or !p1.empty())
    {
         map1 = surf1.shortest_edge_map(p1,-0.1);
         map2 = surf2.shortest_edge_map(p2,-0.1);

         surf1.adjusting_boundary_region(map1.begin() ,map1.end());
         surf2.adjusting_boundary_region(map2.begin() ,map2.end());

         // After adjusting the boundary smoothing is needed, taubin is least volatile.
         surf1.smooth_taubin_region(p1.begin(), p1.end(), 2);
         surf2.smooth_taubin_region(p2.begin(), p2.end(), 2);

         if (iter % 2 == 0)
         {
             p1 = surf2.inside(surf1, p1);
             p2 = surf1.inside(surf2, p2);
         }

         if (iter++ > 60)
             break;
    }
    vertex_vector cp1 = surf1.get_close_points(surf2);
    vertex_vector cp2 = surf2.get_close_points(surf1);

    //seperate_surface as close intersections.
    while (!cp2.empty() or !cp1.empty())
    {
        cp1  = domain.outside(surf1,cp1);
        cp2  = domain.outside(surf2,cp2);

        map1 = surf1.shortest_edge_map(cp1, -0.1); // changed from 0.5
        map2 = surf2.shortest_edge_map(cp2, -0.1);

        surf1.adjusting_boundary_region(map1.begin(), map1.end());
        surf2.adjusting_boundary_region(map2.begin(), map2.end());

        surf1.smooth_taubin_region(cp1.begin(), cp1.end(), 2);
        surf2.smooth_taubin_region(cp2.begin(), cp2.end(), 2);

        cp1 = surf1.get_close_points(surf2);
        cp2 = surf2.get_close_points(surf1);
        if (iter++ > 60)
            break;
    }
}


template< typename CGALSurface>
std::shared_ptr<CGALSurface> morphological_surface_union(CGALSurface& surf1, CGALSurface& surf2)
{
    typedef typename CGALSurface::vertex_vector vertex_vector;
    typedef typename CGALSurface::vertex_descriptor vertex_descriptor;

    std::map< vertex_descriptor, double> map1;
    std::map< vertex_descriptor, double> map2;

    vertex_vector p1 = surf1.get_vertices();
    vertex_vector p2 = surf2.get_vertices();

    p1 = surf2.inside(surf1, p1);
    p2 = surf1.inside(surf2, p2);

    surf1.surface_eval(p1);
    surf2.surface_eval(p2);

    map1 = surf1.shortest_edge_map(p1,1.0);
    map2 = surf2.shortest_edge_map(p2,1.0);

    surf1.adjusting_boundary_region(map1.begin() ,map1.end());
    surf2.adjusting_boundary_region(map2.begin() ,map2.end());

    surf1.smooth_taubin_region(p1.begin(), p1.end(), 2);
    surf2.smooth_taubin_region(p2.begin(), p2.end(), 2);

    std::shared_ptr<CGALSurface> result(new CGALSurface());
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(
            surf1.get_mesh(), surf2.get_mesh(), result->get_mesh());
    return result;
}


class CGALSurface {
    public:
        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
        typedef Kernel::Point_3 Point_3;
        typedef CGAL::Surface_mesh<Point_3> Mesh;
        typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
        typedef Tr::Geom_traits GT;
        typedef Kernel::Vector_3 Vector_3;

        typedef boost::graph_traits< Mesh >::halfedge_descriptor    halfedge_descriptor;
        typedef boost::graph_traits< Mesh >::face_descriptor        face_descriptor;
        typedef boost::graph_traits< Mesh >::vertex_descriptor      vertex_descriptor;
        typedef std::vector< vertex_descriptor >                    vertex_vector;
        typedef CGAL::Side_of_triangle_mesh< Mesh, Kernel > Inside; // declear inside implementation ?
        typedef Mesh::Vertex_index Index;

        typedef CGAL::Polyhedron_3< Kernel > Polyhedron;

        typedef std::vector< Point_3 >  Polyline_3;
        typedef std::vector< Polyline_3 > Polylines;
        typedef Kernel::FT FT;
        typedef std::vector< std::size_t > Face;

        CGALSurface(){} // empty constructor

        CGALSurface(Polyhedron &polyhedron);

        CGALSurface(std::vector<Point_3>& points, std::vector<Face>& faces );

        CGALSurface(const std::string  filename1, const std::string  filename2);

        CGALSurface(const std::string  filename);

        // ---------------------------------------------
        ~CGALSurface(){}

        void split_edges(const double target_edge_length);

        void make_cylinder(const double x0, const double y0, const double  z0, const double x1,
                const double y1, const double z1, const double radius, const int number_of_segments=360) ;

        void make_cone(const double x0, const double y0, const double  z0, const double x1, const double y1,
                const double z1, const double r0, const double r1, const int number_of_segments=360);

        void make_cube(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1);

        void make_sphere(double x0, double y0, double z0, const double r0);

        void clip(const double x1, const double x2, const double x3, const double x4, const double x5,
                const double x6, const bool clip);

        void triangulate_hole();

        Mesh& get_mesh() {return mesh; } // operator overload pybind error due to non-const

        void surface_intersection(CGALSurface& other);

        void surface_difference(CGALSurface& other);

        void surface_union(CGALSurface& other);

        void smooth_taubin(const size_t nb_iter);

        template<typename InputIterator>
        void smooth_taubin_region(InputIterator begin, InputIterator end, const size_t iter);

        void seperate_close_junctures();

        vertex_vector inside(CGALSurface &other, CGALSurface::vertex_vector &points);

        vertex_vector outside(CGALSurface &other, CGALSurface::vertex_vector &points);

        vertex_vector get_vertices();

        std::map< vertex_descriptor, double > seperate_close_surfaces(CGALSurface& other);

        std::shared_ptr< CGALSlice > mesh_slice(double x1, double x2, double x3, double x4);

        std::pair< double, double > span( int direction);

        void surface_eval(vertex_vector &input);

        int collapse_edges(const double stop_ratio);

        int fill_holes();

        bool triangulate_faces();

        void stitch_borders();

        void insert_surface(CGALSurface& surf){mesh += surf.get_mesh(); }

        void clear(){ mesh.clear(); }

        // Global functions
        void isotropic_remeshing(const double target_edge_length,
                const unsigned int nb_iter, const bool protect_border);

        void adjust_boundary(const double c);

        void smooth_laplacian(const double c);

        void fix_close_junctures(const double c);

        std::shared_ptr<CGALSlice> slice(double x1, double x2, double x3, double x4);

        // Local Functions
        void fair(CGALSurface::vertex_vector vector);

        void insert_points(std::vector<Point_3>& points) ;

        // TODO: typedef input itertator ?? 
        template<typename InputIterator>  // -> operation on surface_mesh-> index-based 
        void adjusting_boundary_region(InputIterator begin, InputIterator end, const double c);

        void adjusting_boundary_region(std::map< vertex_descriptor, double >::iterator begin,
                std::map< vertex_descriptor, double >::iterator end);

        template<typename InputIterator>
        void smooth_laplacian_region(InputIterator begin, InputIterator end, const double c);

        vertex_vector get_close_points(CGALSurface &other);

        int num_faces() const {return mesh.number_of_faces(); }

        int num_edges() const { return mesh.number_of_edges(); }

        int num_vertices() const {return mesh.number_of_vertices(); }

        int num_self_intersections();

        bool self_intersections();

        Polylines get_features() { Polylines plines; return plines;} // should be virtual or something ->

        void reconstruct(const double sm_angle = 20.0,
                         const double sm_radius = 100.0,
                         const double sm_distance = 0.25,
                         const double approximation_ratio = 0.02,
                         const double average_spacing_ratio = 5.0);

        template< typename Polyhedron_3>  // Template to address the different kernels in CGAL TODO: return Polyhedron?
        void get_polyhedron(Polyhedron_3 &polyhedron_3) { CGAL::copy_face_graph(mesh, polyhedron_3); }

        void save(const std::string outpath);

        template<typename Implicit_function>
        void implicit_surface(Implicit_function implicit_function,
            const double bounding_sphere_radius,
            const double angular_bound=30.,
            const double radius_bound=0.1,
            const double distance_bound=0.1);

        void write_STL(const std::string filename);

        std::map< vertex_descriptor, double > shortest_edge_map(
                CGALSurface::vertex_vector &vector,double adjustment = -0.5);

    protected:
        Mesh mesh;
};


CGALSurface::CGALSurface(const std::string filename1, const std::string filename2) {
    Mesh mesh1;
    Mesh mesh2;
    utils::load_surface(mesh1, filename1);
    utils::load_surface(mesh2, filename2);
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh1, mesh2, mesh);
}


CGALSurface::CGALSurface(const std::string filename) {
    utils::load_surface(mesh, filename);
}


CGALSurface::CGALSurface(Polyhedron &polyhedron) {
    CGAL::copy_face_graph(polyhedron, mesh);
}


CGALSurface::CGALSurface(std::vector< Point_3 > &points, std::vector< Face > &faces)
{
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, mesh);

    if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
    {
        std::cout<< "reverse_face_orientation"<< std::endl;
        CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
    }
}




//--------------------------------------------
//      CGALSURFACE FUNCTIONS
//--------------------------------------------


void CGALSurface::split_edges(const double target_edge_length) {
    CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
}


void CGALSurface::make_cylinder(const double x0, const double y0, const double z0, const double x1, 
        const double y1, const double z1, const double radius, const int number_of_segments) {
    // TODO :allow tilted cubes
    CGALSurface::make_cone( x0, y0,z0, x1,y1, z1, radius , radius,  number_of_segments);
}


void CGALSurface::make_cone(const double x0, const double y0, const double z0, const double x1,
    const double y1, const double z1, const double r0, const double r1, const int number_of_segments) {
    //TODO: Handle case for r1=0 or r0=0
    Mesh m;

    Index v0 = m.add_vertex(Point_3(x0, y0, z0));
    Index v1 = m.add_vertex(Point_3(x1, y1, z1));

    Index vb;
    Index vt;
    //normalizing vectors ...

    const double n0 = x1 - x0;
    const double n1 = y1 - y0;
    const double n2 = z1 - z0;

    const double l1 = std::sqrt(n0*n0 + n1*n1 + n2*n2);
    const double l2 = std::sqrt((n1 - n2)*(n1 - n2) + (n2 - n0)*(n2 - n0) + (n0 - n1)*(n0 - n1));

    const Vector_3 normal(n0/l1, n1/l1, n2/l1);
    const Vector_3 t1( (n1 - n2)/l2 , (n2 - n0)/l2, (n0 - n1)/l2);
    const Vector_3 t2 = CGAL::cross_product(normal, t1); // normalized ??

    double c  = 360.0 / (double) number_of_segments;        // FIXME: proper casting

    for(int i = 0; i < number_of_segments; ++i) {
        Point_3 pb = Point_3(x0, y0, z0) + t1*r0*std::cos(c*i*CGAL_PI/180) + t2*r0*std::sin(c*i*CGAL_PI/180);
        vb = m.add_vertex(pb);

        Point_3 pt = Point_3(x1, y1, z1) + t1*r1*std::cos(c*i*CGAL_PI/180) + t2*r1*std::sin(c*i*CGAL_PI/180);
        vt = m.add_vertex(pt);

        if (i != 0) {
            m.add_face(v0, vb,Index(vb-2));
            m.add_face(v1, Index(vt-2), vt);

            m.add_face(Index(vt-2), Index(vb-2), vt);
            m.add_face(vb, vt, Index(vb-2));
        }
    }
    m.add_face(Index(0), Index(2), vb);
    m.add_face(Index(1), vt, Index(3));

    m.add_face(vt, vb, Index(3));
    m.add_face(Index(2), Index(3), vb);
    this->mesh = m;
}


void CGALSurface::make_cube(const double x0, const double y0, const double z0, const double x1,
        const double y1, const double z1) {
   // TODO :allow tilted cubes
    assert(x0!=x1);
    assert(y0!=y1);
    assert(z0!=z1);

    Mesh m;
    vertex_descriptor v0 = m.add_vertex(Point_3(x0,y0,z0));
    vertex_descriptor v1 = m.add_vertex(Point_3(x0,y0,z1));
    vertex_descriptor v2 = m.add_vertex(Point_3(x1,y0,z0));

    vertex_descriptor v3 = m.add_vertex(Point_3(x1,y0,z1));
    vertex_descriptor v4 = m.add_vertex(Point_3(x1,y1,z0));
    vertex_descriptor v5 = m.add_vertex(Point_3(x1,y1,z1));

    vertex_descriptor v6 = m.add_vertex(Point_3(x0,y1,z0));
    vertex_descriptor v7 = m.add_vertex(Point_3(x0,y1,z1));

    // Side 1
    m.add_face(v0, v2, v1);
    m.add_face(v3, v1, v2);

    // Side 2
    m.add_face(v3 ,v2, v5);
    m.add_face(v4 ,v5,v2);

    // Side 3
    m.add_face(v4, v6, v5);
    m.add_face(v7, v5, v6);

    // Side 4
    m.add_face(v0, v1, v6);
    m.add_face(v7, v6, v1);

    // Side 5
    m.add_face(v3, v5, v1);
    m.add_face(v7, v1, v5);

    // Side 6
    m.add_face(v0, v6, v2);
    m.add_face(v4, v2, v6);

    this->mesh = m;
}


void CGALSurface::make_sphere(double x0, double y0, double z0, const double r0)
{
    Sphere_wrapper sphere;

    sphere.radius = r0;
    sphere.x0 = x0;
    sphere.y0 = y0;
    sphere.z0 = z0;

    surface_mesher(mesh, sphere.function, x0, y0, z0, r0, 30, r0/5, r0/5);
}


void CGALSurface::clip(const double x1, const double x2, const double x3, const double x4, const double x5,
        const double x6, const bool clip)
{
    // create a plane through point p orthogonal to point v.
    const auto v = Vector_3(x1, x2, x3);    // Surface normal
    const auto p = Point_3(x4, x5, x6);

    CGAL::Polygon_mesh_processing::clip(mesh, Kernel::Plane_3(p, v),
            CGAL::Polygon_mesh_processing::parameters::clip_volume(clip));
}


void CGALSurface::triangulate_hole()
{
    for (const auto h: halfedges(mesh))
    {
        std::vector<face_descriptor> patch_facets;
        if (is_border(h, mesh))
            CGAL::Polygon_mesh_processing::triangulate_hole(mesh, h, back_inserter(patch_facets));
    }
}


void CGALSurface::surface_intersection(CGALSurface &other)
{
    CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh, other.get_mesh(), mesh);
}


void CGALSurface::surface_difference(CGALSurface &other)
{
    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(mesh, other.get_mesh(), mesh);
}


void CGALSurface::surface_union( CGALSurface& other)
{
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh, other.get_mesh(), mesh);
}


void CGALSurface::smooth_taubin(const size_t nb_iter)
{
    for (size_t i = 0; i < nb_iter; ++i)
    {
        this->smooth_laplacian(0.8);
        this->smooth_laplacian(-0.805);
    }
}


int CGALSurface::collapse_edges(const double stop_ratio)
{
    CGAL::Surface_mesh_simplification::Count_ratio_stop_predicate<Mesh> stop(stop_ratio);

    const int r = CGAL::Surface_mesh_simplification::edge_collapse(mesh, stop,
        CGAL::parameters::get_cost(CGAL::Surface_mesh_simplification::Edge_length_cost<Mesh>())
            .get_placement(CGAL::Surface_mesh_simplification::Midpoint_placement<Mesh>()));
    return r;
}


int CGALSurface::fill_holes()
{
    unsigned int nb_holes = 0;
    for(auto h: halfedges(mesh))
    {
        if(is_border(h, mesh))
        {
            std::vector<face_descriptor> patch_facets;
            std::vector<vertex_descriptor> patch_vertices;

            bool success = CGAL::cpp11::get<0>( CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                    mesh, h,
                    std::back_inserter(patch_facets), std::back_inserter(patch_vertices),
                    CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).
                    geom_traits(Kernel())) );

            nb_holes++;
        }
    }
    return nb_holes;
}


bool CGALSurface::triangulate_faces()
{
    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

    for(const auto fit: faces(mesh))
    {
        if (next(next(halfedge(fit, mesh), mesh), mesh) != prev(halfedge(fit, mesh), mesh))
        {
            std::cerr << "Error: non-triangular face left in mesh." << std::endl;
            return false;
        }
    }
    return true;
}


void CGALSurface::stitch_borders()
{
    CGAL::Polygon_mesh_processing::stitch_borders(mesh);

}


void CGALSurface::isotropic_remeshing(const double target_edge_length,
        const unsigned int nb_iter, const bool protect_border)
{
    CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
    CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh), target_edge_length, mesh,
            CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
            .protect_constraints(protect_border));
}


void CGALSurface::fix_close_junctures(const double c)
{
   // name seperate junction and add negtive value as default
   // TODO :  Clean up typedefs, algorithm COULD be faster
   typedef boost::graph_traits<Mesh>::vertex_descriptor                     Point;
   typedef boost::graph_traits<Mesh>::vertices_size_type                    size_type;
   typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type             Vertex_point_pmap;

   // New
   typedef CGAL::Search_traits_3<Kernel>                                    Traits_base;
   typedef CGAL::Search_traits_adapter<Point,Vertex_point_pmap,Traits_base> Traits;
   typedef CGAL::Orthogonal_k_neighbor_search<Traits>                       K_neighbor_search;
   typedef K_neighbor_search::Tree                                          Tree;
   typedef Tree::Splitter                                                   Splitter;
   typedef K_neighbor_search::Distance                                      Distance;

   CGALSurface::vertex_vector results;
   Vertex_point_pmap vppmap = get(CGAL::vertex_point, mesh);

   Tree tree(vertices(mesh).begin(), vertices(mesh).end(), Splitter(), Traits(vppmap));
   Distance tr_dist(vppmap);

   FT distance;
   FT edgeL;
   Point_3 closest;

   bool flag;
   for (const auto vit: mesh.vertices())
   {
        flag = true;

        K_neighbor_search search(tree, mesh.point(vit), 2, 0, true, tr_dist);    //tree.closest_point(point_query); must be second closest
        closest = mesh.point((search.begin() + 1)->first); //closest point is the query point, therefore chose second point 
        Point_3 current = mesh.point(vit);
        distance = CGAL::squared_distance(current, closest);
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(vit), mesh), done(vbegin);
        //ditr_dist.inverse_of_transformed_distance(search.begin()->second);
        do
        {
            edgeL = CGAL::squared_distance(current, mesh.point(*vbegin));
            if (distance >= edgeL)
            {
                flag=false;
                break;
            }
            *vbegin++;
        }
        while(vbegin != done);

        if (flag)
            results.push_back(vit);
   }

   adjusting_boundary_region(results.begin(), results.end(), c);
   //fair(results);
}


std::shared_ptr<CGALSlice> CGALSurface::slice(double x1, double x2, double x3, double x4)
{
    typedef Kernel::Point_2 Point_2;
    typedef std::vector<Point_2> Polyline_2;
    typedef std::vector<Polyline_2> Polylines_2;
    //-------------------------------------------------------------
    // Intersection Polylines
    //-------------------------------------------------------------
    CGAL::Polygon_mesh_slicer<Mesh, Kernel> slicer(mesh);
    Polylines polylines_3D;

    slicer(Kernel::Plane_3(x1, x2, x3, x4), std::back_inserter(polylines_3D));

    Kernel::Plane_3 plane = Kernel::Plane_3(x1, x2, x3, x4);
    Vector_3 n =  plane.orthogonal_vector();

    Vector_3 e1 = plane.base1();
    Vector_3 e2 = plane.base2();

    //-------------------------------------------------------------
    // Polylines to 2D
    //-------------------------------------------------------------

    std::vector<std::vector<Point_2>> polylines_2;
    for (const auto pol: polylines_3D)
    {
        std::vector<Point_2> result;
        std::vector<Point_2> polyline_2;
        for (const auto pit: pol)
        {
            polyline_2.push_back(plane.to_2d(pit));
        }
        polylines_2.push_back(polyline_2);
    }

    std::shared_ptr<CGALSlice> slice(new CGALSlice(polylines_2)); 
    slice->set_plane(plane); 
    return slice;
}


void CGALSurface::fair(CGALSurface::vertex_vector svector)
{
    CGAL::Polygon_mesh_processing::fair(mesh, svector);
}


void CGALSurface::insert_points(std::vector<Point_3>& points)
{
    // template iterator ?? ADD POINT ineffective, but could be implemented
    for(const auto it: points)
        mesh.add_vertex(it);
}


void CGALSurface::adjust_boundary(const double c)
{
    Mesh::Vertex_range::iterator vb = mesh.vertices().begin(), ve = mesh.vertices().end();
    CGALSurface::adjusting_boundary_region(vb, ve, c);
}


void CGALSurface::smooth_laplacian(const double c)
{
    Mesh::Vertex_range::iterator vb = mesh.vertices().begin(), ve = mesh.vertices().end();
    CGALSurface::smooth_laplacian_region(vb, ve, c);
}


template<typename InputIterator>
void CGALSurface::adjusting_boundary_region(InputIterator begin , InputIterator  end, const double c)
{
    /* typedef typename CGAL::Vertex_around_target_circulator<Mesh> HV_const_circulator; */
    std::vector<std::pair<vertex_descriptor, Point_3> > smoothed;
    for ( ; begin != end; ++begin)
    {
        Vector_3 delta= CGAL::Polygon_mesh_processing::compute_vertex_normal(*begin,mesh);
        Point_3 p = mesh.point(*begin) + c*delta;
        smoothed.push_back(std::make_pair(*begin, p));
    }

    for (const auto s: smoothed)
        mesh.point(s.first) = s.second;
}


template<typename InputIterator >
void CGALSurface::smooth_laplacian_region(InputIterator begin, InputIterator end, const double c)
{
    std::vector<std::pair<vertex_descriptor, Point_3> > smoothed;

    for ( ; begin != end; ++begin)
    {
        Point_3 current = mesh.point(*begin);
        Vector_3 delta=CGAL::NULL_VECTOR;
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(*begin),mesh), done(vbegin);

        do
        {
            delta += Vector_3(mesh.point(*vbegin) - current);
            *vbegin++;
        }
        while(vbegin != done);
        auto p = current + c*delta/mesh.degree(*begin);
        smoothed.push_back(std::make_pair(*begin, p));
    }

    for (const auto s: smoothed)
        mesh.point(s.first) = s.second;
}


int CGALSurface::num_self_intersections()
{
    std::vector< std::pair<face_descriptor, face_descriptor> > intersected_tris;
    CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(intersected_tris));
    return intersected_tris.size();  // Could actually return the triangles themselves
}


bool CGALSurface::self_intersections()
{
    return CGAL::Polygon_mesh_processing::does_self_intersect(mesh,
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));
}


void CGALSurface::reconstruct(const double sm_angle, const double sm_radius, const double sm_distance,
                              const double approximation_ratio, const double average_spacing_ratio)
{
    reconstruct::poisson_reconstruction(mesh, sm_angle, sm_radius, sm_distance,
            approximation_ratio, average_spacing_ratio);
}


void CGALSurface::save(const std::string outpath)
{
    std::ofstream out(outpath);
    out << mesh;
    out.close();
}


template<typename Implicit_function>
void CGALSurface::implicit_surface(Implicit_function implicit_function,
            const double bounding_sphere_radius,
            const double angular_bound,
            const double radius_bound,
            const double distance_bound)
{
    surface_mesher(mesh, implicit_function, bounding_sphere_radius, angular_bound,
            radius_bound, distance_bound);
}


std::map< CGALSurface::vertex_descriptor, double > CGALSurface::shortest_edge_map(
        CGALSurface::vertex_vector &vertex_vector, double adjustment)
{
    std::map< vertex_descriptor, double > results;
    for (const auto v_it: vertex_vector)
    {
        Point_3 current = mesh.point(v_it);
        CGAL::Vertex_around_target_circulator< Mesh > vbegin(mesh.halfedge(v_it), mesh), done(vbegin);

        FT min_edge = FT(100);
        do
        {
            FT temp = CGAL::squared_distance(current, mesh.point(*vbegin));
            if (temp < min_edge )
            {
                min_edge = temp;
            }
         *vbegin++;
        } while (vbegin != done);

        results[v_it] = adjustment*static_cast< double >(CGAL::sqrt(min_edge));
    }
    return results;
}


CGALSurface::vertex_vector CGALSurface::get_vertices()
{
    // CLEANER
    vertex_vector result;
    for (const auto v_it: mesh.vertices())
    {
        result.emplace_back(v_it);
    }
    return result;
}


CGALSurface::vertex_vector CGALSurface::inside(CGALSurface &other, CGALSurface::vertex_vector &points)
{
    vertex_vector result;
    CGALSurface::Inside inside_poly2(mesh);

    for (const auto v_it: points)
    {
        CGAL::Bounded_side res = inside_poly2(other.get_mesh().point(v_it));
        if (res == CGAL::ON_BOUNDED_SIDE or res == CGAL::ON_BOUNDARY)
            result.emplace_back(v_it);
    }
    return result;
}


CGALSurface::vertex_vector CGALSurface::outside(CGALSurface &other, CGALSurface::vertex_vector &points)
{
    vertex_vector result;
    CGALSurface::Inside inside_poly2(mesh);
    for (const auto v_it: points)
    {
        CGAL::Bounded_side res = inside_poly2(other.get_mesh().point(v_it));
        if (res == CGAL::ON_UNBOUNDED_SIDE)
            result.emplace_back(v_it);
    }
    return result;
}


void CGALSurface::seperate_close_junctures() // TODO: based on distance to surface
{
    typedef boost::graph_traits<Mesh>::vertex_descriptor                     Point;
    typedef boost::graph_traits<Mesh>::vertices_size_type                    size_type;
    typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type             Vertex_point_pmap;

    typedef CGAL::Search_traits_3<Kernel>                                    Traits_base;
    typedef CGAL::Search_traits_adapter<Point,Vertex_point_pmap,Traits_base> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>                       K_neighbor_search;
    typedef K_neighbor_search::Tree                                          Tree;
    typedef Tree::Splitter                                                   Splitter;
    typedef K_neighbor_search::Distance                                      Distance;

    std::map<vertex_descriptor,double> results;
    Vertex_point_pmap vppmap = get(CGAL::vertex_point, mesh);

    Tree tree(vertices(mesh).begin(), vertices(mesh).end(), Splitter(), Traits(vppmap));

    Distance tr_dist(vppmap);

    FT distance;
    FT edgeL;
    Point_3 closest;
    bool flag;
    for (vertex_descriptor v_it : mesh.vertices())
    {
        flag = true;
         //tree.closest_point(point_query); must be second closest
        K_neighbor_search search(tree, mesh.point(v_it), 2, 0, true, tr_dist);

        // The closest point is the query point, therefore choose second point as closest.
        closest = mesh.point((search.begin() + 1)->first);

        Point_3 current = mesh.point(v_it);
        distance = CGAL::squared_distance(current, closest);
        CGAL::Vertex_around_target_circulator< Mesh > vbegin(mesh.halfedge(v_it), mesh), done(vbegin);
        do
        {
            edgeL = CGAL::squared_distance(current, mesh.point(*vbegin));
            if (distance >= edgeL)
            {
                flag = false;
                break;
            }
        *vbegin++;
        } while (vbegin != done);

        if (flag)
        {
            // both vertices are affected changes from 0.5
            results[v_it] = -0.2*static_cast<double>(CGAL::sqrt(distance));
        }
    }

    adjusting_boundary_region(results.begin(), results.end());
    smooth_taubin(2);
}


CGALSurface::vertex_vector CGALSurface::get_close_points(CGALSurface &other)
{
    typedef boost::graph_traits<Mesh>::vertex_descriptor                     Point;
    typedef boost::graph_traits<Mesh>::vertices_size_type                    size_type;
    typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type             Vertex_point_pmap;

    typedef CGAL::Search_traits_3<Kernel>                                    Traits_base;
    typedef CGAL::Search_traits_adapter<Point,Vertex_point_pmap,Traits_base> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>                       K_neighbor_search;
    typedef K_neighbor_search::Tree                                          Tree;
    typedef Tree::Splitter                                                   Splitter;
    typedef K_neighbor_search::Distance                                      Distance;

    CGALSurface::vertex_vector results;
    Vertex_point_pmap vppmap = get(CGAL::vertex_point,other.get_mesh());

    Tree tree(vertices(other.get_mesh()).begin(), vertices(other.get_mesh()).end(), Splitter(), Traits(vppmap));
    Distance tr_dist(vppmap);

    FT distance;
    FT edgeL;
    Point_3 closest;
    bool flag;
    for (vertex_descriptor v_it : mesh.vertices())
    {
        flag = true;
        // every point as a close point
        K_neighbor_search search(tree, mesh.point(v_it), 2, 0, true, tr_dist);

        closest = other.get_mesh().point((search.begin())->first);

        Point_3 current = mesh.point(v_it);
        distance = CGAL::squared_distance(current, closest );
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(v_it),mesh), done(vbegin);
        do
        {
            edgeL = CGAL::squared_distance(current, mesh.point(*vbegin));
            if (distance >= edgeL)
            {
                flag=false;
                break;
            }
            *vbegin++;
        } while (vbegin != done);

        if (flag)
            results.push_back(v_it); // both vertices are affected
    }
    return results;
}


std::map< CGALSurface::vertex_descriptor, double > CGALSurface::seperate_close_surfaces(CGALSurface& other)
{
    // TODO: based on distance to surface
    typedef boost::graph_traits<Mesh>::vertex_descriptor                     Point;
    typedef boost::graph_traits<Mesh>::vertices_size_type                    size_type;
    typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type             Vertex_point_pmap;

    typedef CGAL::Search_traits_3<Kernel>                                    Traits_base;
    typedef CGAL::Search_traits_adapter<Point,Vertex_point_pmap,Traits_base> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>                       K_neighbor_search;
    typedef K_neighbor_search::Tree                                          Tree;
    typedef Tree::Splitter                                                   Splitter;
    typedef K_neighbor_search::Distance                                      Distance;

    std::map< vertex_descriptor, double > results;
    Vertex_point_pmap vppmap = get(CGAL::vertex_point,other.get_mesh());

    Tree tree(vertices(other.get_mesh()).begin(), vertices(other.get_mesh()).end(), Splitter(), Traits(vppmap));
    Distance tr_dist(vppmap);

    FT distance;
    FT edgeL;
    Point_3 closest;
    bool flag;
    for (vertex_descriptor v_it: mesh.vertices())
    {
        flag = true;
        //tree.closest_point(point_query); must be second closest
        K_neighbor_search search(tree, mesh.point(v_it), 2,0,true,tr_dist);

        // The closest point is the query point, therefore choose second point as closest.
        closest = mesh.point((search.begin()+1)->first);

        Point_3 current = mesh.point(v_it);
        distance = CGAL::squared_distance(current, closest);
        CGAL::Vertex_around_target_circulator< Mesh > vbegin(mesh.halfedge(v_it), mesh), done(vbegin);
        do
        {
            edgeL = CGAL::squared_distance(current, mesh.point(*vbegin));
            if (distance >= edgeL)
            {
                flag = false;
                break;
            }
            *vbegin++;
        } while (vbegin != done);

        if (flag)
            results[v_it] = -0.2*static_cast<double>(CGAL::sqrt(distance));  // both vertices are affected
    }

    return results;
}


template< typename InputIterator >
void CGALSurface::smooth_taubin_region(InputIterator begin, InputIterator end, const size_t nb_iter)
{
    for (size_t i = 0; i < nb_iter; ++i) {
        this->smooth_laplacian_region(begin, end, 0.8);
        this->smooth_laplacian_region(begin, end, -0.805);
    }
}


void CGALSurface::surface_eval(CGALSurface::vertex_vector &input)
{
    typedef typename CGAL::Vertex_around_target_circulator<Mesh> HV_const_circulator;
    int size = input.size();

    Vector_3 direction;
    for (int i = 0; i < size; ++i)
    {
        Vector_3 n1 = CGAL::Polygon_mesh_processing::compute_vertex_normal(input[i], mesh);
        direction += n1;
        CGAL::Vertex_around_target_circulator< Mesh > vbegin(mesh.halfedge(input[i]), mesh), done(vbegin);
        do
        {
            Vector_3 n2 = CGAL::Polygon_mesh_processing::compute_vertex_normal(*vbegin, mesh);

            if (n1*n2 > CGAL::sqrt(n1.squared_length())*CGAL::sqrt(n2.squared_length())*0.50)
            {
                vertex_descriptor tar = *vbegin;

                if(std::find(input.begin(), input.end(), tar) == input.end())
                {
                    input.push_back(tar);
                    ++size;
                }
            }

            *vbegin++;
        } while (vbegin != done);
    }
}


#endif
