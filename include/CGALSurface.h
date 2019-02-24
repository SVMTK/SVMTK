#ifndef CGALSurface_H


#define CGALSurface_H

#ifndef BOOST_PARAMETER_MAX_ARITY
# define BOOST_PARAMETER_MAX_ARITY 12
#endif
 


// Local
#include "read_polygons_STL.h"
#include "surface_mesher.h"
#include "reconstruct_surface.h"




// BOOST
#include <boost/foreach.hpp>
#include <CGAL/boost/graph/copy_face_graph.h>

//  STL
#include <assert.h>  
#include <iterator>
#include <vector>
#include <CGAL/IO/STL_reader.h>






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


#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Polyline_simplification_2/simplify.h>

//---------------------------------------------------------
// CGAL- Polygon_mesh_processing
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
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>

#include <CGAL/Side_of_triangle_mesh.h> 



//--------------------------------------------
//      UTILITY
//--------------------------------------------
template< typename Mesh> // seperate because -> typedef mesh in CGALSurface
bool load_surface(const std::string file, Mesh& mesh) // new name load polygons if load fails
{
   
  std::ifstream input(file);

  std::string extension = file.substr(file.find_last_of(".")+1);//TODO: FIX boost linking problem
  if (!input)
  {
    std::cerr << "Cannot open file " << std::endl;
    return false;
  }

  typedef typename Mesh::Point Point_3;
  
  
  std::vector<Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;

  if ( extension=="off")
  {
     std::cout<< "reading off" << std::endl;
     if (!CGAL::read_OFF(input, points, polygons))
     {
         std::cerr << "Error parsing the OFF file " << std::endl;
         return false;
     }
     std::cout<< "finished" << std::endl;
  }
  else if ( extension=="stl")
  {
     if (!read_polygons_STL(input, points, polygons)) 
     {
         std::cerr << "Error parsing the STL file " << std::endl;
         return false;
     }
  }
  else
  {
     std::cerr << "Error unkown file" << std::endl;
     return false;
  }

  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);


  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,mesh);

  if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
  {
   std::cout<< "reverse_face_orientation"<< std::endl;

    CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
  }
  return true;
}



template< typename CGALSurface>  // FIXME
void surface_overlapp(CGALSurface& surf1 , CGALSurface& surf2 , double c,    int max_iter  = 300)
{
    typedef typename CGALSurface::vertex_vector vertex_vector;

    vertex_vector surf1points;
    vertex_vector surf2points;

    surf1points = surf1.points_inside(surf2);
    surf2points = surf2.points_inside(surf1);

    int iter =0;

    while (  !surf1points.empty() and  !surf2points.empty() )
    {
            surf1.fair(surf1points);
            surf2.fair(surf2points);

            surf1points= surf1.points_inside(surf2 , surf1points);
            surf2points = surf2.points_inside(surf1, surf2points);

            surf1.adjusting_boundary_region(surf1points.begin() ,surf1points.end(), c);
            surf2.adjusting_boundary_region(surf2points.begin() ,surf2points.end(), c);

            //std::cout<< surf1points.size() + surf2points.size() << std::endl;
            iter++;
            if( iter >max_iter)
                 break;
              
    }

}




class CGALSurface
{
   public:
  
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel; // Change to exact -> copy face_ graph to inexact ??
    typedef Kernel::Point_3 Point_3;
    typedef CGAL::Surface_mesh<Point_3> Mesh;
    typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
    typedef Tr::Geom_traits GT;
    typedef Kernel::Vector_3 Vector_3;

    typedef boost::graph_traits<Mesh>::halfedge_descriptor    halfedge_descriptor;
    typedef boost::graph_traits<Mesh>::face_descriptor        face_descriptor;
    typedef boost::graph_traits<Mesh>::vertex_descriptor      vertex_descriptor; // 
    typedef std::vector<vertex_descriptor>                    vertex_vector;
    typedef CGAL::Side_of_triangle_mesh<Mesh,Kernel> Inside; // declear inside implementation ? 
    typedef Mesh::Vertex_index Index;

    typedef CGAL::Polyhedron_3<Kernel> Polyhedron; // just for 

    typedef std::pair<Point_3, Vector_3> RPwn;
    typedef std::vector<Point_3>  Polyline_3; 
    typedef std::vector<Polyline_3> Polylines; 
    typedef Kernel::FT FT;

    CGALSurface(){} // empty constructor

    CGALSurface(const std::string  filename1, const std::string  filename2);
    CGALSurface(const std::string  filename);

    CGALSurface(Polyhedron &polyhedron); 
    // TODO: remove 1 one these

    // ---------------------------------------------
    ~CGALSurface(){}

        void split_edges(double  target_edge_length);

        void make_cylinder( double x0, double y0, double  z0,  double x1, double y1, double z1, double radius,  int number_of_segments=360) ;

        void make_cone( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0 , double r1,  int number_of_segments=360) ;

        void make_cube( double x0, double y0, double  z0,  double x1, double y1, double z1);

        void make_sphere( double x0, double y0, double  z0, double r0);

        void clip(double x1,double x2, double x3 ,double x4, bool clip);
        
        void triangulate_hole();

        Mesh& get_mesh() {return mesh;} // operator overload pybind error due to non-const
    
        void surface_intersection( CGALSurface& other);

        void surface_difference( CGALSurface& other);

        void surface_union( CGALSurface& other);

        void smooth_taubin(const size_t nb_iter); 

        int collapse_edges(const double stop_ratio);

        int fill_holes();

        bool triangulate_faces();
      
        void stitch_borders();
         
        bool insert_surface(CGALSurface& surf){mesh+=surf.get_mesh();} // 

        void clear(){ mesh.clear();}
        // Global functions 
        void isotropic_remeshing(double target_edge_length, unsigned int nb_iter, bool protect_border);
    
        void adjust_boundary(const double c);

        void smooth_laplacian(const double c, int iter);

        void fix_close_junctures(double c);

        void mesh_slice(double x1,double x2, double x3 ,double x4, std::string outpath) ;

        // Local Functions
        void fair(CGALSurface::vertex_vector vector); 

        void insert_points(std::vector<Point_3>& points) ;

        // TODO: typedef input itertator ?? 
        template<typename InputIterator>  // -> operation on surface_mesh-> index-based 
        void adjusting_boundary_region(InputIterator begin , InputIterator  end, const double c);

        template<typename InputIterator>
        void smooth_laplacian_region(InputIterator begin , InputIterator end ,const double c);

        //-----------------------------------------------------------------------------------
        //TODO: Add query and overload ? ?
        vertex_vector points_inside(CGALSurface& other); // Query overlaod 

        vertex_vector points_outside(CGALSurface& other);

        vertex_vector  points_inside(CGALSurface& other,CGALSurface::vertex_vector &points);

        void reconstruct_surface(const double sm_angle,const double sm_radius,const double sm_distance) ;

        int num_faces() const {return mesh.number_of_faces();}

        int num_edges() const { return mesh.number_of_edges();}

        int num_vertices() const {return mesh.number_of_vertices();}

        int num_self_intersections();
   
        bool self_intersections(); 


        Polylines get_features() { Polylines plines; return plines;} // should be virtual or something ->
        
        void reconstruct( double sm_angle = 20.0,
                          double sm_radius = 100.0,
                          double sm_distance = 0.25,
                          double approximation_ratio = 0.02,
                          double average_spacing_ratio = 5.0);

        template< typename Polyhedron_3>  // Template to address the different kernels in CGAL TODO: return Polyhedron?
        void get_polyhedron(Polyhedron_3 &polyhedron_3 ){CGAL::copy_face_graph(mesh,polyhedron_3);}

        void save(const std::string outpath);

        template<typename Implicit_function>  
        void implicit_surface(Implicit_function implicit_function,
             double bounding_sphere_radius,
             double angular_bound=30.,
             double radius_bound=0.1 ,
             double distance_bound=0.1   );
        //void preprocess(double target_edge_length, int nb_iter); // TODO:REMOVE or implement


        //void surface_extension(double x0 , double y0 , double z0); provide point -> closest point on surface computes surface normal and creates cylinder with point as midpoint
        void write_STL(const std::string filename);


   //private :
   protected:
       Mesh mesh;


};


void CGALSurface::write_STL(const std::string filename)
{

    std::ofstream file(filename);
    file.precision(6);
    Vector_3 n;

    file << "solid "<< filename << std::endl;

    for (auto f : mesh.faces() )
    {

       n = CGAL::Polygon_mesh_processing::compute_face_normal(f,mesh);


       file << "facet normal " << n.x() << " " << n.y()   << " " << n.z() <<std::endl;
       file << "outer loop"<< std::endl;
       CGAL::Vertex_around_face_iterator<Mesh> vbegin, vend;
       for(boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(f), mesh); vbegin != vend; ++vbegin)
       {
          file << "\t" << "vertex " << mesh.point(*vbegin).x() <<" "<< mesh.point(*vbegin).y() <<" "<< mesh.point(*vbegin).z()  << std::endl;  
            
       }
       file <<"endloop" << std::endl;
       file <<"endfacet"<< std::endl;
    }
    file << "endsolid"  << std::endl;

}




template<typename Implicit_function>
void CGALSurface::implicit_surface(Implicit_function implicit_function,
             double bounding_sphere_radius,
             double angular_bound,
             double radius_bound,
             double distance_bound)
{
     surface_mesher(mesh,implicit_function,bounding_sphere_radius,angular_bound,radius_bound, distance_bound);
}


template<typename InputIterator >
void CGALSurface::adjusting_boundary_region(InputIterator begin , InputIterator  end, const double c)
{
  typedef typename CGAL::Vertex_around_target_circulator<Mesh> HV_const_circulator;
  std::vector<std::pair<vertex_descriptor, Point_3> > smoothed; 

  for ( ; begin != end; ++begin)
  {
      Vector_3 delta= CGAL::Polygon_mesh_processing::compute_vertex_normal(*begin,mesh); 
      Point_3 p = mesh.point(*begin) + c*delta;
      smoothed.push_back(std::make_pair(*begin, p));
  }
  for (std::pair<vertex_descriptor, Point_3> s : smoothed)
  {
      mesh.point(s.first) = s.second;
  }
}
template<typename InputIterator >
void CGALSurface::smooth_laplacian_region(InputIterator begin , InputIterator end ,const double c)
{
  typedef typename CGAL::Vertex_around_target_circulator<Mesh> HV_const_circulator;
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
      }while(vbegin!=done);

      Point_3 p = current + c*delta/mesh.degree(*begin);  // WHY C
      smoothed.push_back(std::make_pair(*begin, p));
  }
  for (std::pair<vertex_descriptor, Point_3> s : smoothed)
  {
      mesh.point(s.first) = s.second;
  }
}




//--------------------------------------------
//      CGALSURFACE FUNCTIONS
//--------------------------------------------


CGALSurface::CGALSurface(const std::string  filename1, const std::string  filename2)
{ 
    Mesh mesh1;
    Mesh mesh2; 
    load_surface(filename1,mesh1);
    load_surface(filename2,mesh2);
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh1, mesh2,mesh);

}
CGALSurface::CGALSurface(const std::string filename)
{
    load_surface(filename,mesh);
}

// Boolean operation FIXME to overloaded operators ?
void CGALSurface::surface_intersection( CGALSurface& other)
{
     CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh, other.get_mesh(),mesh) ;
}
void CGALSurface::surface_difference( CGALSurface& other)
{
     CGAL::Polygon_mesh_processing::corefine_and_compute_difference(mesh , other.get_mesh(),mesh);
}

void CGALSurface::surface_union( CGALSurface& other)
{
     CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh, other.get_mesh(),mesh);
}


void CGALSurface::smooth_taubin(const size_t nb_iter) {
    for (size_t i = 0; i < nb_iter; ++i) {
        this->smooth_laplacian(0.8,1);
        this->smooth_laplacian(-0.805,1);
    }
}

int CGALSurface::num_self_intersections() {
    std::vector< std::pair<face_descriptor, face_descriptor> > intersected_tris;
    CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(intersected_tris));
    return intersected_tris.size();  // Could actually return the triangles themselves -> return facets
}

int CGALSurface::collapse_edges(const double stop_ratio) {
    
    CGAL::Surface_mesh_simplification::Count_ratio_stop_predicate<Mesh> stop(stop_ratio);

    const int r = CGAL::Surface_mesh_simplification::edge_collapse(
        mesh,
        stop,
        CGAL::parameters::get_cost(CGAL::Surface_mesh_simplification::Edge_length_cost<Mesh>())
            .get_placement(CGAL::Surface_mesh_simplification::Midpoint_placement<Mesh>()));

    return r;
}

void CGALSurface::mesh_slice(double x1,double x2, double x3 ,double x4, std::string outpath)
{

     // NOTE :: NOT WORKING
     typedef Kernel::Point_2 Point_2;
     typedef std::vector<Point_2> Polyline_2;
     typedef std::vector<Polyline_2> Polylines_2;

     


     // Intersection Polylines
     //-------------------------------------------------------------
     CGAL::Polygon_mesh_slicer<Mesh, Kernel> slicer(mesh); 
     Polylines polylines_3D;

     slicer(Kernel::Plane_3(x1, x2, x3, x4), std::back_inserter(polylines_3D));

     //-------------------------------------------------------------
     // Basis in slice 
     //-------------------------------------------------------------
     Vector_3 n(x1,x2,x3); 
     
     n=n/CGAL::sqrt(n.squared_length());

     Vector_3 b1( 1,0,0);
     
     Vector_3 e1 = b1 - (b1*n)*n;    

     e1=e1/CGAL::sqrt(e1.squared_length());
    
     Vector_3 e2 = CGAL::cross_product(n,e1);

     //-------------------------------------------------------------
     // Polylines in 2D
     //-------------------------------------------------------------
  
     std::vector<std::vector<Point_2>> polylines_2;
     
     for ( auto  pol = polylines_3D.begin(); pol != polylines_3D.end(); ++pol ) 
     {
         std::vector<Point_2> points;
         for ( auto pit = pol->begin(); pit != pol->end(); ++pit)
         {
               points.push_back(Point_2( pit->x()*e1.x() + pit->y()*e1.y() + pit->z()*e1.z() , pit->x()*e2.x() + pit->y()*e2.y() + pit->z()*e2.z() ));
         }
         polylines_2.push_back(points);
     }
   
     // Polylines in 2D -> simplify polylines in 2D



     typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb_2;
     typedef CGAL::Delaunay_mesh_face_base_2<Kernel> Fb_2;
     typedef CGAL::Triangulation_data_structure_2<Vb_2, Fb_2> Tds_2;
     typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds_2,CGAL::Exact_predicates_tag> CDT_2;
     typedef CGAL::Delaunay_mesh_size_criteria_2<CDT_2> Criteria;
     typedef CGAL::Delaunay_mesher_2<CDT_2, Criteria> Mesher;

     typedef CGAL::Delaunay_mesh_size_criteria_2<CDT_2> Criteria;
     typedef CGAL::Delaunay_mesher_2<CDT_2, Criteria> Mesher;

     CDT_2 cdt;


 
     for ( auto pol = polylines_2.begin(); pol != polylines_2.end(); ++pol)
     {
          cdt.insert_constraint(pol->begin(), pol->end() );
          break;
     }

     //CGAL::Polyline_simplification_2::simplify(ct, Cost(), Stop(0.8));




     //CGAL::refine_Delaunay_mesh_2(ct, Criteria());

    
     Mesher mesher(cdt);
     mesher.refine_mesh();

     for(CDT_2::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
     {
         if (!fit->is_in_domain())
         {
            cdt.delete_face(fit);
         }
     } 

     std::ofstream out(outpath);
     CGAL::export_triangulation_2_to_off(out,cdt);
  


     
}




bool CGALSurface::self_intersections()
{
    return CGAL::Polygon_mesh_processing::does_self_intersect(mesh,
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));
}

void CGALSurface::isotropic_remeshing(double target_edge_length, unsigned int nb_iter, bool protect_border)
{

     CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
     CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh),
                              target_edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
                              .protect_constraints(protect_border));

}

void CGALSurface::clip(double x1,double x2, double x3 ,double x4, bool clip)
{
   CGAL::Polygon_mesh_processing::clip(mesh,  Kernel::Plane_3(x1,x2,x3,x4), CGAL::Polygon_mesh_processing::parameters::clip_volume(clip));

}

void CGALSurface::triangulate_hole()
{

    
    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh))
    {
        std::vector<face_descriptor>  patch_facets;
        if(is_border(h,mesh))
        {
           CGAL::Polygon_mesh_processing::triangulate_hole(mesh, h,back_inserter(patch_facets));
        }
    }

}
int CGALSurface::fill_holes()
{
    unsigned int nb_holes = 0;
    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh))
    {
      if(is_border(h,mesh))
      {
        std::vector<face_descriptor>  patch_facets;
        std::vector<vertex_descriptor> patch_vertices;
        bool success = CGAL::cpp11::get<0>(
          CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                    mesh,
                    h,
                    std::back_inserter(patch_facets),
                    std::back_inserter(patch_vertices),
          CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).
                    geom_traits(Kernel())) );

        std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
        std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
        std::cout << "  Is fairing successful: " << success << std::endl;
        nb_holes++;
      }
    }

    std::cout << std::endl;
    std::cout << nb_holes << " holes have been filled" << std::endl;

    return nb_holes;
}
bool CGALSurface::triangulate_faces()
{
    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

    BOOST_FOREACH(face_descriptor fit, faces(mesh))
      if (next(next(halfedge(fit, mesh), mesh), mesh)
          !=   prev(halfedge(fit, mesh), mesh))
        std::cerr << "Error: non-triangular face left in mesh." << std::endl;
    return true;
}
void CGALSurface::stitch_borders()
{
    CGAL::Polygon_mesh_processing::stitch_borders(mesh);

}

CGALSurface::vertex_vector CGALSurface::points_inside(CGALSurface& other)
{
  vertex_vector result;
  CGALSurface::Inside inside_poly2(other.get_mesh());

  for ( vertex_descriptor v_it : mesh.vertices())
  {
      CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
      if (res == CGAL::ON_BOUNDED_SIDE or res == CGAL::ON_BOUNDARY){result.push_back(v_it);}
  }
  return result;
}

CGALSurface::vertex_vector CGALSurface::points_inside(CGALSurface& other,CGALSurface::vertex_vector &points)
{
  vertex_vector result;
  CGALSurface::Inside inside_poly2(other.get_mesh());

  for ( vertex_descriptor v_it : points )
  {
      CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
      if (res == CGAL::ON_BOUNDED_SIDE or res == CGAL::ON_BOUNDARY){result.push_back(v_it);}
  }
  points.clear();
  return result;
}



CGALSurface::vertex_vector CGALSurface::points_outside(CGALSurface& other)
{
  vertex_vector result;
  CGALSurface::Inside inside_poly2(other.get_mesh());

  for ( vertex_descriptor v_it : mesh.vertices())
  {
      CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
      if (res == CGAL::ON_UNBOUNDED_SIDE){result.push_back(v_it);}
  }
  return result;
}




void CGALSurface::adjust_boundary(const double c)
{
    Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    CGALSurface::adjusting_boundary_region(vb, ve,c);

}
void CGALSurface::smooth_laplacian(const double c, int iter)
{
    Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    for ( int i = 0 ; i< iter ; ++i)
    {
        CGALSurface::smooth_laplacian_region(vb, ve,c);
    }
}

void CGALSurface::make_cube( double x0, double y0, double  z0,  double x1, double y1, double z1) // TODO :allow tilted cubes
{

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

void CGALSurface::make_cylinder( double x0, double y0, double  z0,  double x1, double y1, double z1, double radius,  int number_of_segments) // TODO :allow tilted cubes
{
      
        CGALSurface::make_cone( x0, y0,z0, x1,y1, z1, radius , radius,  number_of_segments ) ;

}
CGALSurface::CGALSurface(Polyhedron &polyhedron) {
    CGAL::copy_face_graph(polyhedron, mesh);
}


void CGALSurface::make_cone( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0, double r1,  int number_of_segments) // TODO :allow tilted cubes
{
        //TODO: Handle case for r1=0 or r0=0
        Mesh m;

    
        Index v0 = m.add_vertex(Point_3(x0,y0,z0));
        Index v1 = m.add_vertex(Point_3(x1,y1,z1));
        
        Index vb;
        Index vt;
        //normalizing vectors ...

        double n0 = x1-x0;
        double n1 = y1-y0;
        double n2 = z1-z0;

        double l1 = std::sqrt(n0*n0+n1*n1+n2*n2);

        Vector_3 normal( n0/l1 , n1/l1, n2/l1 );

        double l2 = std::sqrt( (n1-n2)*(n1-n2) + (n2-n0)*(n2-n0) + ( n0 -n1)*(n0-n1) );

        Vector_3 t1( (n1-n2)/l2 , (n2-n0)/l2, (n0-n1)/l2 ); 

        Vector_3 t2 = CGAL::cross_product(normal,t1); // normalized ?? 

        double c  = 360.0/(double)number_of_segments;
     

        for(int i = 0; i < number_of_segments; ++i) 
        {
              
           Point_3 pb= Point_3(x0,y0,z0)  + t1*r0*std::cos(c*i*CGAL_PI/180) + t2*r0*std::sin(c*i*CGAL_PI/180); 

           vb = m.add_vertex(pb);

           Point_3 pt = Point_3(x1,y1,z1)  + t1*r1*std::cos(c*i*CGAL_PI/180) + t2*r1*std::sin(c*i*CGAL_PI/180); 

           vt = m.add_vertex(pt);

           if ( i!=0)
           {  

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





struct sphere_wrapper{
        
     public:
        static double radius;
        static double x0;
        static double y0;
        static double z0;
        
        static double function(double x, double y , double z) { return  (x-x0)*(x -x0) +  (y-y0)*(y -y0) + (z-z0)*(z -z0) -radius*radius; }

        

};

// TODO : Sphere uses a implicit function, and can be extended to other implict fucntion with smooth boundaries.  
double sphere_wrapper::radius = 0;
double sphere_wrapper::x0 = 0;
double sphere_wrapper::y0 = 0;
double sphere_wrapper::z0 = 0;

void CGALSurface::make_sphere( double x0, double y0, double  z0,double r0 ) 
{

  sphere_wrapper sphere;
   
  sphere.radius = r0;
  sphere.x0=x0;
  sphere.y0=y0;
  sphere.z0=z0;

  surface_mesher(mesh,sphere.function,x0,y0,z0,r0,30,r0/5,r0/5);   

}

void CGALSurface::split_edges(double  target_edge_length)
{
     CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
}




void CGALSurface::save(const std::string outpath)
{    
     std::string extension = outpath.substr(outpath.find_last_of(".")+1);
     std::ofstream out(outpath);
     if ( extension=="off")
     {
        out << mesh;
        out.close();
     } 
     else if ( extension=="stl")
     {
       write_STL(outpath);
     }


}

void CGALSurface::fair(CGALSurface::vertex_vector vector)
{
     CGAL::Polygon_mesh_processing::fair(mesh,vector);
}





void CGALSurface::insert_points(std::vector<Point_3>& points)  // template iterator ?? ADD POINT ineffective, but could be implemented
{
     for( std::vector<Point_3>::iterator it=points.begin(); it!= points.end(); ++it){mesh.add_vertex(*it);}
}

void CGALSurface::reconstruct( double sm_angle,
                               double sm_radius,
                               double sm_distance,
                               double approximation_ratio,
                               double average_spacing_ratio )
{ 

poisson_reconstruction(mesh,sm_angle,
                            sm_radius,
                            sm_distance,
                            approximation_ratio,
                            average_spacing_ratio);

}


void CGALSurface::fix_close_junctures(double c) // name seperate junction and add negtive value as default
{
   // TODO :  Clean up typedefs, algorithm COULD be faster
   typedef boost::graph_traits<Mesh>::vertex_descriptor                     Point;
   typedef boost::graph_traits<Mesh>::vertices_size_type                    size_type;
   typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type             Vertex_point_pmap;

   // new 
   typedef CGAL::Search_traits_3<Kernel>                                    Traits_base;
   typedef CGAL::Search_traits_adapter<Point,Vertex_point_pmap,Traits_base> Traits;
   typedef CGAL::Orthogonal_k_neighbor_search<Traits>                       K_neighbor_search;
   typedef K_neighbor_search::Tree                                          Tree;
   typedef Tree::Splitter                                                   Splitter;
   typedef K_neighbor_search::Distance                                      Distance;

 
   CGALSurface::vertex_vector results;
   Vertex_point_pmap vppmap = get(CGAL::vertex_point,mesh);

   Tree tree(vertices(mesh).begin(),
            vertices(mesh).end(),
            Splitter(),
            Traits(vppmap)
   );

   Distance tr_dist(vppmap);

   FT distance;
   FT edgeL;
   Point_3 closest;
   bool flag;
   for (vertex_descriptor v_it : mesh.vertices())
   {
        flag = true;

        K_neighbor_search search(tree, mesh.point(v_it), 2,0,true,tr_dist); //tree.closest_point(point_query); must be second closest
        
        closest = mesh.point((search.begin()+1)->first); //closest point is the query point, therefore chose second point 
        
        Point_3 current = mesh.point(v_it);
        distance = CGAL::squared_distance(current, closest );
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(v_it),mesh), done(vbegin);
        //ditr_dist.inverse_of_transformed_distance(search.begin()->second);
        do
        {
  
           edgeL = CGAL::squared_distance(current, mesh.point(*vbegin) ); 
   
           if ( distance >= edgeL ) 
           {
              flag=false;
            break;
           }
         *vbegin++;
        }while(vbegin!=done);

        if (flag){results.push_back(v_it);}

   }
   
   adjusting_boundary_region(results.begin(),results.end(),c);
   
   //fair(results); 
   
}


#endif
