// Copyright (C) 2018-2021 Lars Magnus Valnes and Jakob Schreiner
//
// This file is part of Surface Volume Meshing Toolkit (SVM-TK).
//
// SVM-Tk is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SVM-Tk is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SVM-Tk.  If not, see <http://www.gnu.org/licenses/>.

#ifndef Surface_H
#define Surface_H


#ifndef BOOST_PARAMETER_MAX_ARITY
# define BOOST_PARAMETER_MAX_ARITY 12
#endif

#define _default_prratio 1.73
//#define _prratio 2.0 //1.414

/* --Includes -- */
#include "Errors.h"  

/* -- boost-- */
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <CGAL/Surface_mesh.h>

/* -- CGAL 3D Convex Hulls  -- */
#include <CGAL/convex_hull_3.h>

/* -- CGAL Principal Component Analysis -- */
#include <CGAL/centroid.h>

/* -- CGAL dD Spatial Searching  -- */
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>


/* -- CGAL Triangulated Surface Mesh Skeletonization -- */
#include <CGAL/extract_mean_curvature_flow_skeleton.h>


/* -- CGAL IO -- */
#include <CGAL/IO/STL.h> 
#include <CGAL/IO/OFF.h> 

/* -- CGAL Polygon Mesh Processing -- */
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h> 
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path_traits.h>
#include <CGAL/Surface_mesh_shortest_path.h>
//#include <CGAL/Surface_mesh_default_triangulation_3.h>
/* -- CGAL Surface mesh simplification -- */
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <CGAL/boost/graph/Euler_operations.h>

//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>

/* -- CGAL AABB -- */
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

/* -- CGAL Triangulated Surface Mesh Segmentation -- */

#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Kernel/global_functions.h>


#include <CGAL/Polyhedron_3.h>

// typedef double (Function) (const Point_3&);
// typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function*> Labeling_function;
/**
 * \class 
 * Cost function used to collapse edges with 
 * zero cost of removal, i.e. edges in a straight line. 
 */
template<class TM_>
class Cost_stop_predicate
{
  public:
    typedef TM_ TM;

    Cost_stop_predicate(double threshold) : thres(threshold) {}

    template< typename F, typename Profile>
    bool operator()( F const & aCurrentCost, Profile const & profile, std::size_t aic,std::size_t acc) const
    {
      return static_cast<double>(aCurrentCost) > thres;
    }
  private :
    double thres;
};


template<typename Point_3>
Point_3 vector_centroid( std::vector< Point_3> points ) 
{
       return CGAL::centroid( points.begin(), points.end());
}
    
/**
 * @note polyline implies some structure 
 */    
template<typename Point_3> 
void ordering_polyline_3(std::vector< Point_3> &points, std::string ordering) 
{
      if ( ordering == "x-axis" ) 
      {
            if ( points.front().x() > points.back().x() )
                 std::reverse( points.begin() , points.end() ) ;
      }
      else if ( ordering == "y-axis" ) 
      {
            if ( points.front().y() > points.back().y() )
                 std::reverse( points.begin() , points.end() ) ;
      }
      else if ( ordering == "z-axis" ) 
      {
            if ( points.front().z() > points.back().z() )
                 std::reverse( points.begin() , points.end() ) ;
      }
      else 
        return;      
}



template<typename Point_3> 
std::vector<std::vector<Point_3>> ordering_polylines_3( std::vector<std::vector<Point_3>> &polylines, std::string ordering) 
{
      // Make more efficient
      std::vector<std::pair<double,int>> min_values; 
      int iter = 0;   
      for( auto polyline : polylines ) 
      {
           ordering_polyline_3( polyline,ordering);
           if      ( ordering == "x-axis" ) 
                min_values.push_back( std::make_pair<double,int>( CGAL::to_double( polyline.front().x()), iter++) );
           else if ( ordering == "y-axis" ) 
                min_values.push_back( std::make_pair<double,int>( CGAL::to_double( polyline.front().y()), iter++) );
           else if ( ordering == "z-axis" ) 
                min_values.push_back( std::make_pair<double,int>( CGAL::to_double( polyline.front().z()), iter++) );
           else 
               continue;
      }
      
      if (min_values.empty())
          return polylines;
          
      struct sort_by_first_element {
                  bool operator()(const std::pair<double,int> & a, const std::pair<double,int> & b)
                       { return( a.first < b.first );}           
      };
   
      std::sort(min_values.begin(), min_values.end(), sort_by_first_element() ); 

      std::vector<std::vector<Point_3>> output;      
    
      for ( auto srt : min_values ) 
      {
            ordering_polyline_3(polylines[srt.second],ordering);
            output.push_back(polylines[srt.second]);
      } 
      return output;               
}

template< typename InputIterator> 
double length_polyline_3( InputIterator begin , InputIterator end)
{
  double length = 0.0;
  for(; begin != end; ++begin)
     length += CGAL::to_double( CGAL::sqrt(CGAL::squared_distance(*begin, *(begin+1))));        
  return length; 
} 

template< typename Point_3 > 
double distance_to_polyline( std::vector<Point_3> points  , Point_3 query)
{
  double min_length = CGAL::to_double( CGAL::sqrt(CGAL::squared_distance(points[0], query)));
  double length;
  for( auto point : points )
  {
        length = CGAL::to_double( CGAL::sqrt(CGAL::squared_distance(point, query)));    
        if ( min_length > length)
              min_length = length;
  }
  return min_length; 
} 



// Euclidean
template<typename Point_3> 
std::vector<std::vector<Point_3>> get_longest_polyline( std::vector<std::vector<Point_3>> &polylines)
{
      std::vector<Point_3> output;      
      double max_length = 0;
      double length;
      for( auto polyline : polylines ) 
      {
          length = length_polyline_3(polyline.begin(),polyline.end());
          if ( max_length < length ) 
          {
              max_length = length;
               output = polyline;
          }
      }
      return output ;
} 

template<typename Point_3> 
std::vector<std::vector<Point_3>> get_polyline_near_point( std::vector<std::vector<Point_3>> &polylines, Point_3 point)
{
      std::vector<Point_3> output;      
      double max_length = 0;
      double length;
      for( auto polyline : polylines ) 
      {
          length = distance_to_polyline(polyline.begin(),polyline.end(), point );
          if ( max_length < length ) 
          {
              max_length = length;
               output = polyline;
          }
      }
      return output ;
} 

// DocString: Surface
/**
 * \class Surface 
 *
 * The SVMTK Surface class is used to create, store and manipulate triangulated surfaces in 3D. 
 * For auxilary purpose, most of the CGAL declarations are done inside the Surface class.
 *
 * CGAL is implemented with different kernels, that have different properties, and this 
 * class is implemnted with the Exact predicates inexact constructions kernel.
 * @see(CGAL::Exact_predicates_inexact_constructions_kernel)[https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html] 
 * 
 * SVMTK Surface class is used to handle operations related to triangualted 
 * surfaces in 3D CGAL Surface Mesh .
 * @see (CGAL::Surface_mesh<Point_3>)[https://doc.cgal.org/latest/Surface_mesh/index.html]
 *   
 * 
 * It should be noted that most operations can also be used with CGAL
 * Polyhedral Surface.
 * @see(CGAL::Polyhedron_3<Kernel>)[https://doc.cgal.org/latest/Polyhedron/index.html]
 *
 */ 
class Surface
{
  public:
   typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel; 
   
   typedef Kernel::Point_3 Point_3;
   typedef Kernel::Plane_3 Plane_3;
   typedef Kernel::FT FT; 
   typedef Kernel::Ray_3 Ray;  
   typedef Kernel::Segment_3 Segment_3; 
   typedef Kernel::Line_3 Line_3;         
   typedef Kernel::Sphere_3 Sphere_3; 
   typedef Kernel::Vector_3 Vector_3;
   typedef Kernel::Triangle_3 Triangle_3;
   typedef Kernel::Point_2 Point_2;
    
   typedef CGAL::Surface_mesh<Point_3> Mesh;
   typedef CGAL::Side_of_triangle_mesh<Mesh,Kernel> Inside; 
   typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;   
   typedef CGAL::Vertex_around_target_circulator<Mesh> HV_const_circulator; 
   //typedef CGAL::Face_around_target_circulator<Mesh> HF_const_circulator; 

   typedef boost::graph_traits<Mesh>::edge_descriptor             edge_descriptor;
   typedef boost::graph_traits<Mesh>::halfedge_descriptor         halfedge_descriptor;
   typedef boost::graph_traits<Mesh>::face_descriptor             face_descriptor;
   typedef boost::graph_traits<Mesh>::vertex_descriptor           vertex_descriptor; 
   typedef boost::graph_traits<Mesh>::vertices_size_type          size_type;
   typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type   Vertex_point_pmap; 
  
   typedef std::vector<Point_3>               Polyline;     // Polylines is a point vector with a specific sequence of points.
   typedef std::vector<Point_3>               point_vector;   
   typedef std::vector<Polyline>              Polylines;
   typedef std::vector<Vector_3>              vertex_disp;             
   
   typedef std::vector<std::size_t>           Face;
   
   typedef std::vector<vertex_descriptor>     vertex_vector;   
   typedef std::set<vertex_descriptor>        vertex_set;         

   typedef std::vector<face_descriptor>       face_vector; 
   typedef std::set<face_descriptor>          face_set;            
   
   typedef Mesh::Vertex_index                 Index;

   typedef std::map<vertex_descriptor,Vector_3> vertex_vector_map;
   typedef std::map<vertex_descriptor,double> vertex_scalar_map;     
       
   typedef CGAL::Search_traits_3<Kernel>                                                Traits_base;
   typedef CGAL::Search_traits_adapter<vertex_descriptor,Vertex_point_pmap,Traits_base> Traits;
   typedef CGAL::Orthogonal_k_neighbor_search<Traits>                                   K_neighbor_search;
   typedef K_neighbor_search::Tree                                                      Tree;
   typedef Tree::Splitter                                                               Splitter;
   typedef K_neighbor_search::Distance                                                  Distance;

   typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
   typedef CGAL::AABB_traits_3<Kernel, Primitive> AABB_Traits_3;
   typedef CGAL::AABB_tree<AABB_Traits_3> AABB_Tree;
       
   typedef std::optional<AABB_Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
   
   typedef boost::property_map<Mesh,CGAL::edge_is_feature_t>::type EIFMap; 
   typedef boost::property_map<Mesh, CGAL::face_patch_id_t<int> >::type PIMap;


   typedef CGAL::Mean_curvature_flow_skeletonization<Mesh> Skeletonization;
   typedef Skeletonization::Skeleton                             Skeleton;
   typedef Skeleton::vertex_descriptor                           Skeleton_vertex;

   typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Mesh> shortest_path_traits;
   typedef CGAL::Surface_mesh_shortest_path<shortest_path_traits> Surface_mesh_shortest_path;
   typedef Surface_mesh_shortest_path::Face_location Face_location;
   
   
   // DocString: Surface
   /**
    * @brief Constructs an empty SVMTK Surface object.
    */
   Surface() : _stratio(5.0e-3) ,  _prratio(_default_prratio), _smreduc(0.8) {} 

   /**
    * @brief Takes a CGAL polyhedron surface and copies the structure over to CGAL surface mesh.
    * @param polyhedron  defined CGAL::Polyhedron_3<Kernel> Polyhedron; 
    */
    Surface(Polyhedron_3 &polyhedron) :   _stratio(5.0e-3) ,  _prratio(_default_prratio), _smreduc(0.8)
    {
      CGAL::copy_face_graph(polyhedron, mesh);
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);          
    }

   // DocString: Surface
   /**
    * @brief Constructs a SVMTK Surface object with surface from file. 
    *        Current fileformats: off, stl  
    * @param filename the string path to surface to load.
    */
    Surface(const std::string filename) :    _stratio(5.0e-3) ,  _prratio(_default_prratio), _smreduc(0.8)
    {
      load_surface(filename);
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);   


    }

   // DocString: Surface
   /**
    * @brief Constructs a copy of a SVMTK Surface object.
    */
    Surface(const std::shared_ptr<Surface> surf) :   _stratio(5.0e-3) ,  _prratio(_default_prratio), _smreduc(0.8)
    {
       this->mesh=surf.get()->get_mesh(); 
       CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);  

    } 
    
   // DocString: Surface
   /**
    * @brief Constructs a copy of a SVMTK Surface object.
    */
    Surface(Mesh mesh) :    _stratio(5.0e-3) ,  _prratio(_default_prratio), _smreduc(0.8)
    {
       this->mesh=mesh; 
       CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);    
    } 
   
    
   /**
    * @brief Constructs a SVMTK Surface object with surface using points and connections between points.
    * @param points coordinates of vertices  
    * @param faces connections of  vertices
    */
    Surface(std::vector<Point_3>& points, std::vector<Face>& faces) : _stratio(5.0e-3) ,  _prratio(_default_prratio), _smreduc(0.8)
    {
       CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces);
       CGAL::Polygon_mesh_processing::repair_polygon_soup(points, faces);
       CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces,mesh);
       if(CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
       {
         CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
       }
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);       

    }

    Surface(std::vector<Point_3>& points, std::vector<std::array<int,3>>& faces) :    _stratio(5.0e-3) ,  _prratio(_default_prratio), _smreduc(0.8)
    {
       CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces);
       CGAL::Polygon_mesh_processing::repair_polygon_soup(points, faces);
       CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces,mesh);
       if(CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
       {
         CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
       }
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);       

    }

   ~Surface(){}
    
  /**
   * @brief Copy Operator
   */  
   Surface &operator=(Surface &other) { this->mesh = other.mesh; return *this; }

  /**
   * @brief Loads triangulated surfaces with exentsion off or stl.
   *
   * Loads surface mesh as a polygon soup, and 
   * reassembles it according to CGAL structure
   * of faces and vertices.
   * 
   * @param file string of the fileaneme.
   * @returns true if successful.
   */
   bool load_surface(const std::string file)
   {
     std::vector<Point_3> points;
     std::vector< std::vector<std::size_t> > polygons;
     std::ifstream input(file);
     std::string extension = file.substr(file.find_last_of(".") + 1);
     if( !input )  
     {
       throw InvalidArgumentError("Cannot open file");
       return false;
     }
     if( extension=="off" )
     {
       if( !CGAL::IO::read_OFF(input, points, polygons) )
       {
         throw InvalidArgumentError("Error parsing the OFF file."); 
         return false;
       }
     }
     else if( extension=="stl" )
     {
       if( !CGAL::IO::read_STL(input, points, polygons) )
       {
         throw InvalidArgumentError("Error parsing the STL file.");
         return false;
       }
     }
     else
     {
       throw InvalidArgumentError("Error unknown file.");
       return false;
     }
     CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
     CGAL::Polygon_mesh_processing::repair_polygon_soup(points, polygons);
     CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,this->mesh);

     if( CGAL::is_closed(this->mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(this->mesh)) )
         CGAL::Polygon_mesh_processing::reverse_face_orientations(this->mesh);
     return true;
   }
  

  /** 
   *
   */
   void add_points(std::vector<Point_3> points)
   {
     for( Point_3 point : points)
         mesh.add_vertex(point);
   }


  /** 
   *
   */  
   void add_point(Point_3 point)
   {
      mesh.add_vertex(point);
   }
   
  // DocString: get_close_vertex_displacement     
  // 
  /**
   * Experimental: test
   *   Needs improvment on performance and convergence issue.
   * @brief Computes the displacement for each vertices parallel to the shortest vector to vertex  
   *        related to another surface.
   *  Implicit vertex displacement : 
   *               p_f =  p_0 + adj*normal + beta*Sum^N (p_i - p_f )  | + beta*N*p_0 - beta*N*p_0
   * ( 1 + N*beta) p_f = ( 1 + beta*N)p_0 + adj*normal + beta*Sum^N ( p_i - p_0)  
   * p_f - p_0 = (  adj*normal + beta*Sum^N ( p_i - p_0) )/( 1 + beta*N) 
   * 
   * Assumption
   * As long as beta > 0, then displacment p_f - p_0 --> 0 .  
   * @param double multiplier for the displacement.
   * @param double smoothing factor for the displacement.      
   * @returns returns the centeroid point of the surface
   * F
   */    
   template<int A=0> 
   vertex_vector_map get_close_vertex_displacement(Surface other,vertex_vector &mvertices, const double adjustment, const double smoothing)
   {
      Vertex_point_pmap vppmap = get(CGAL::vertex_point,other.get_mesh());  

      Tree tree(vertices(other.get_mesh()).begin(),
                vertices(other.get_mesh()).end(),
                Splitter(),
                Traits(vppmap));
 
       Distance tr_dist(vppmap);
       FT edge_length, min_edge;
       Point_3 closest;
       Vector_3 direction,normal, edge ;
       
       std::map< vertex_descriptor, Vector_3> vertex_displacement;

       double avg_edge_length = average_edge_length();

       double N, distance;
       Vector_3 disp;
       for(auto vit = mvertices.begin(); vit!= mvertices.end(); vit++)
       {
          K_neighbor_search search(tree, mesh.point(*vit), 2,0,true,tr_dist); 
          
          closest = other.get_mesh().point((search.begin()+A)->first);
          
          Point_3 current = mesh.point(*vit);
          
          normal =  CGAL::Polygon_mesh_processing::compute_vertex_normal(*vit,mesh);

          if ( current == closest)
               direction=normal;
          else
          { 
            // Vector(a, b) -> b-a
            direction =  Vector_3(current,closest); 
            distance  = CGAL::sqrt( direction.squared_length() );  
            direction/=distance;
          }

          Vector_3 delta=CGAL::NULL_VECTOR;
          CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(*vit), mesh), done(vbegin);
          N=0;
          min_edge = FT(avg_edge_length);
          do
          {
             N++;
             edge   = Vector_3(current,mesh.point(*vbegin));
             delta += edge;
             edge_length = CGAL::sqrt(edge.squared_length());
             if( edge_length < min_edge )
                 min_edge = edge_length; 
                 
             *vbegin++;
          }while(vbegin!=done);

          disp =  (-adjustment*min_edge*direction + smoothing*delta)/(1 + abs(smoothing*N));  
               
          vertex_displacement[*vit] = disp;
      }   
      return vertex_displacement;
   }

  // DocString: get_vertex_displacement     
  /**
   * @brief Computes the displacement for vertices, which is parallel to the surface normal direction.
   *        
   * 
   * @param double multiplier for the displacement.
   * @param double smoothing factor for the displacement.
   * @returns returns the centeroid point of the surface
   */  
   vertex_vector_map get_vertex_displacement(vertex_vector &mvertices, const double adjustment, const double smoothing)
   {
      std::map< vertex_descriptor, Vector_3> vertex_displacement;

      double avg_edge_length = average_edge_length();
      double N; 
      Vector_3 normal, edge, disp ;
      FT edge_length, min_edge;

      for(auto vit = mvertices.begin(); vit!= mvertices.end(); vit++)
      {
          Point_3 current = mesh.point(*vit);
          normal =  CGAL::Polygon_mesh_processing::compute_vertex_normal(*vit,mesh);
          Vector_3 delta=CGAL::NULL_VECTOR;
          CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(*vit), mesh), done(vbegin);
          N=0;
          min_edge = FT(avg_edge_length);
          do
          {
             edge   = Vector_3(current,mesh.point(*vbegin));
             delta += edge;
             edge_length = CGAL::sqrt(edge.squared_length());
             if( edge_length < min_edge )
                 min_edge = edge_length; 
             *vbegin++;N++;
          }while(vbegin!=done);
          
          disp =  (adjustment*min_edge*normal + smoothing*delta)/(1+abs(smoothing*N)); 
          vertex_displacement[*vit] = disp;
      }
      return vertex_displacement;
   }  

  // DocString: get_dir_vertex_displacement     
  /**
   * @brief Computes the displacement for vertices, which is parallel to the surface normal direction.
   *        
   * 
   * @param double multiplier for the displacement.
   * @param double smoothing factor for the displacement.
   * @returns returns the centeroid point of the surface
   */  
   vertex_vector_map get_constant_vertex_displacement(vertex_vector &mvertices, Vector_3 disp, const double smoothing)
   {
      FT edge_length, min_edge;
      Vector_3 direction,normal, edge ;
      std::map< vertex_descriptor, Vector_3> vertex_displacement;
       
      double N;
      double avg_edge_length = average_edge_length();
      
      for(auto vit = mvertices.begin(); vit!= mvertices.end(); vit++)
      {
          Point_3 current = mesh.point(*vit);
          normal =  CGAL::Polygon_mesh_processing::compute_vertex_normal(*vit,mesh);

          Vector_3 delta=CGAL::NULL_VECTOR;
          CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(*vit), mesh), done(vbegin);
          N=0;
          min_edge = FT(avg_edge_length);
          do
          {
             N++;
             edge   = Vector_3(current,mesh.point(*vbegin));
             delta += edge;
             edge_length = CGAL::sqrt(edge.squared_length());
             if( edge_length < min_edge )
                  min_edge = edge_length; 
             *vbegin++;
          }while(vbegin!=done);

          vertex_displacement[*vit] = ( min_edge*disp + smoothing*delta)/(1+abs(smoothing*N));
      }
      return vertex_displacement;
   }  

  // DocString: add_adjacent_vertices     
  /**
   * @brief Sets the vertex displacement, and adjusts adjacent vertices with similar surface normal.
   *
   * @param vertex_displacement; map with vertices and displacements
   * @param threshold; max angle for similar surface normal.
   * @returns returns the centeroid point of the surface
   */ 

   void add_adjacent_vertices(vertex_vector &verticies, double threshold = 40.00)    
   {
      Vector_3 n1,n2;
      double cos_angle = std::cos(threshold*CGAL_PI/180.0);
      vertex_set nvert; 
      for( auto vit : verticies) 
      {
          n1 =  CGAL::Polygon_mesh_processing::compute_vertex_normal(vit,mesh);
          CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(vit), mesh), done(vbegin);
          do
          {
            n2 =  CGAL::Polygon_mesh_processing::compute_vertex_normal(*vbegin,mesh);         
            if( std::find(verticies.begin(), verticies.end(),*vbegin ) == verticies.end() and  n1*n2>cos_angle ) 
                nvert.insert(*vbegin);
           *vbegin++;
         }while(vbegin!=done);
      }  
      std::copy(nvert.begin(), nvert.end(), std::back_inserter(verticies));
   }
   
   /**
    * @brief Sets the vertex displacement.
    *
    * @param vertex_displacement; map with vertices and displacements
    * @returns none
    */
   void set_vertices(vertex_vector_map &vertex_displacement)    
   {
      for( auto vit : vertex_displacement) 
      {
          mesh.point( vit.first) += vit.second;
      }
   } 
   
   void set_vertex( vertex_descriptor vd, Point_3 point)
   {
        mesh.point( vd) = point;
   }
   
   
  // DocString: set_adjacent_vertices     
  /**
   * @brief Sets the vertex displacement, and adjusts adjacent vertices with similar surface normal.
   *
   * @param vertex_displacement; map with vertices and displacements
   * @param threshold; max angle for similar surface normal.
   * @returns none
   */ 
   void set_adjacent_vertices(vertex_vector_map &vertex_displacement, double threshold = 40.00)    
   {
      Vector_3 n1,n2;
      double cos_angle = std::cos(threshold*CGAL_PI/180.0);
      for( auto vit : vertex_displacement) 
      {
          mesh.point( vit.first) += vit.second;
          n1 =  CGAL::Polygon_mesh_processing::compute_vertex_normal(vit.first,mesh);
     
          CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(vit.first), mesh), done(vbegin);
          do
          {
            n2 =  CGAL::Polygon_mesh_processing::compute_vertex_normal(*vbegin,mesh); 
                    
            if( vertex_displacement.find(vit.first) == vertex_displacement.end() and n1*n2>cos_angle ) 
            {
                mesh.point(*vbegin)+= 0.5*Vector_3(mesh.point(*vbegin),  mesh.point( vit.first) ); // # TODO test with different weight
            }
           *vbegin++;
         }while(vbegin!=done);
      }  
   }
  
  // DocString: centeroid     
  /**
   * @brief Computes the centeroid point of the surface
   * @returns returns the centeroid point of the surface
   */    
   Point_3 centeroid()
   {   
      return CGAL::Polygon_mesh_processing::centroid(mesh);
   }      
   
  // DocString: volume    
  /**
   * @brief Computes the volume enclosed by the surface.
   * @returns the volume enclosed by the surface.
   */    
   double volume()
   {  
      if( !does_bound_a_volume() ) 
          throw PreconditionError("Surface does not bound volume.");
      return CGAL::to_double(CGAL::Polygon_mesh_processing::volume(mesh));
   }
  
  // DocString: area
  /**
   * @brief Computes the area of the surface.
   * @returns the area of the surface
   */
   double area()
   {
      return CGAL::to_double(CGAL::Polygon_mesh_processing::area(mesh));
   } 
  
  /**
   * @brief Throws error if the class variable mesh is empty
   * @returns true if mesh is non empty, otherwise false 
   */
   bool assert_non_empty_mesh()
   { 
      if(mesh.is_empty())
         throw EmptyMeshError("Surface is empty.");
      return true;
   }  
  
  // DocString: average_edge_length
  /**
   * @brief Computes the average edge length in the stored mesh object. 
   * @returns the average edge length
   */
   double average_edge_length()
   {
      double sum = 0; 
      for(edge_descriptor e : mesh.edges())
      {
         halfedge_descriptor he = mesh.halfedge(e);
         sum+=CGAL::to_double( CGAL::Polygon_mesh_processing::edge_length(he,mesh));     
      }
      return sum/mesh.number_of_edges();
   } 
  
  /**
   * @brief Returns the mesh.
   * @returns the mesh of the Surface object.
   */        
   Mesh& get_mesh() 
   {   
      return mesh;
   }
      
  // DocString: clear   
  /**
   * @brief Clears mesh
   */    
   void clear()
   { 
      mesh.clear();
   }
 
   // DocString: num_faces
  /**
   * @brief Returns the number of faces in the surface. 
   * @returns the number of faces in the surface
   */       
   int num_faces() const 
   {
      return mesh.number_of_faces();
   }
 
  // DocString: num_edges
  /**
   * @brief Returns the number of edges in the surface
   * @returns the number of edges in the surface
   */     
   int num_edges() const 
   {
       return mesh.number_of_edges();
   }   
 
   // DocString: num_vertices 
  /**
   * @brief Returns the number of vertices in the surface
   * @returns the number of vertices in the surface
   */     
   int num_vertices() const 
   {
      return mesh.number_of_vertices();
   }

  /**
   * @brief Copies the CGAL Surface Mesh to CGAL Polyhedron_3
   * @param polyhedron_3 CGAL Polyhedron_3 object that is copies the stored mesh.
   */      
   template< typename Polyhedron_3>  
   void get_polyhedron(Polyhedron_3 &polyhedron_3)
   { 
       assert_non_empty_mesh(); 
       CGAL::copy_face_graph(mesh,polyhedron_3);
   }    

  // DocString: implicit_surface    
  /** 
   * @brief  Creates a surface mesh based on an implicit function
   * 
   * @see [Surface_mesher](https://doc.cgal.org/latest/Surface_mesher/index.html)
   * @tparam
   * @tparam an implicit function that takes Cartesian coordinates.
   *          The function has a boundary defined as  
   *         f(x,y,z)=0  and the interior defined as f(x,y,z) < 0 
   * @param implitict function  
   * @param bounding sphere that encloses the mesh construction. , and the
   * @param angular_bound for the minimum facet angle in degrees
   * @param radius_bound  for the radius of the surface Delaunay balls
   * @param distance_bound for center-center distances.
   */

   template<typename Domain, typename Implicit_function>
   void implicit_surface(Implicit_function implicit_function, double mesh_resolution, double bounding_radius, double error_bound)
   {     
   
     auto bounding = Kernel::Sphere_3(CGAL::ORIGIN, FT(bounding_radius*bounding_radius));   
     std::shared_ptr<Domain>   domain( new Domain(implicit_function, bounding,  error_bound));
     domain->create_surface_mesh(mesh_resolution, bounding_radius);
     auto surface = domain->template get_boundary<Surface>(1);
     this->mesh = surface->get_mesh();

    }
     
     
     
     
         





  /** 
   * @brief Takes a vector of vertices and finds cluster of adjacent vertices with normals determined 
   * by  
   *       $ normal*n_2  > cos_angle $ 
   * with normal = Avg(n_i)   
   *
   * @note Smoothing causes most adjacent vertex normals to be consider within the cluster. Therefore, 
   * it is more robust to take the average of all vertex normal corresponding to the input.
   *  
   * @param[in,out] vertices a vector of vertices.  
   * @param angle_in_degree threshold cosinus angle between two vectors 
   */
   void get_normal_vector_cluster(vertex_vector &vertices, double angle_in_degree=36.87)  
   {
      //assert_non_empty_mesh();
      //CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh); 
      double cos_angle = std::cos(angle_in_degree*CGAL_PI/180.0);
      int size = vertices.size();
      Vector_3 normal; 
      for(auto vit : vertices) 
         normal = normal + CGAL::Polygon_mesh_processing::compute_vertex_normal(vit,mesh); 
         
      normal = normal.direction().vector(); 
      for(int i=0; i<size; ++i)
      {
        HV_const_circulator vbegin(mesh.halfedge(vertices[i]),mesh), done(vbegin);
        do
        {
          Vector_3 n2 = CGAL::Polygon_mesh_processing::compute_vertex_normal(*vbegin,mesh);               
          if( normal*n2 > cos_angle )
          {    
            if(std::find(vertices.begin(), vertices.end(), *vbegin) == vertices.end()) 
            {
              vertices.push_back(*vbegin); 
              ++size;
            } 
          } 
          *vbegin++;
        }while(vbegin!=done);
      }
   }
   
  // DocString: reconstruct   
  /**
   * @brief Reconstruct surface mesh. 
   * 
   * Reconstruct a surface based on a CGAL surface mesh object with points using CGAL poisson_reconstruction algorithm.
   * @see [poisson_reconstruction](https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/index.html)
   * @param bounding_sphere_radius indicating the radius of a sphere bounding the meshing operations
   * @param angular_bound bounds for the minimum facet angle in degrees.
   * @param radius_bound bound for the minimum for the radius of the surface Delaunay balls and the center-center distances respectively.
   * @param distance_bound bound for the minimum center-center distances respectively.
   * @overload  
   */
   void reconstruct(double angular_bound=20, double radius_bound=0.1, double distance_bound=0.1)
   { 
      assert_non_empty_mesh();
      //poisson_reconstruction(*this, angular_bound, radius_bound, distance_bound);
   }

   void reconstruct(std::string filename, double angular_bound=20, double radius_bound=0.1, double distance_bound=0.1)
   { 

      //poisson_reconstruction(*this, filename, angular_bound, radius_bound, distance_bound);
      assert_non_empty_mesh(); 
   }
       
  /**
   * @brief Segments the surface, called Domain::boundary_segmentation function.
   *
   * Uses CGAL function sharp edges segmentation that marks sharp edges and segments facet restricted by marked edges.   
   * @see [sharp_edges_segmentation](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__detect__features__grp.html) 
   *
   * @param angle_in_degree sharp edges defined as cosinus angle between two vectors exceed angle_in_degree.
   * @param nb_of_patch_plus_one used as the initial value to mark the surface segmentations. 
   * @returns vector of facets points represented as triangles and associated tags. 
   */
   //std::vector<std::pair<Triangle_3 , std::pair<int,int>>> surface_segmentation(int nb_of_patch_plus_one=1, double angle_in_degree=85)
   std::map<std::pair<int,int>, std::vector<Triangle_3>>  surface_segmentation(int nb_of_patch_plus_one=1, double angle_in_degree=85)
   {
     EIFMap eif = get(CGAL::edge_is_feature, mesh);
    
     Mesh::Property_map<face_descriptor, std::pair<int,int> > patch_id_map;
     Mesh::Property_map<vertex_descriptor,std::set<std::pair<int,int> > > vertex_incident_patch_map;                      

     patch_id_map = mesh.add_property_map<face_descriptor,std::pair<int, int> >("f:pid",std::pair<int,int>()).first; 
     vertex_incident_patch_map = mesh.add_property_map<vertex_descriptor,std::set<std::pair<int, int> > >("f:vip",std::set<std::pair<int, int> >()).first;
        
     CGAL::Polygon_mesh_processing::sharp_edges_segmentation(mesh, angle_in_degree, eif,patch_id_map,
                                     CGAL::Polygon_mesh_processing::parameters::first_index(nb_of_patch_plus_one)
                                    .vertex_incident_patches_map(vertex_incident_patch_map));

     std::vector<std::pair<Triangle_3,std::pair<int,int>>> Tri2tagvec;
     
     std::map<std::pair<int,int>, std::vector<Triangle_3>> ndaf;
   
     Vertex_point_pmap vpm = get(CGAL::vertex_point,mesh);

     for(face_descriptor f : mesh.faces())
     {
         halfedge_descriptor he = mesh.halfedge(f);
         Point_3 p1 = get(vpm,mesh.source(he));
         he = mesh.next(he);
         Point_3 p2 = get(vpm,mesh.source(he));
         he = mesh.next(he);
         Point_3 p3 = get(vpm,mesh.source(he));
         Triangle_3 tri(p1,p2,p3); 
         Tri2tagvec.push_back(std::make_pair(tri,get(patch_id_map,f)));    
          
         if ( ndaf.find(get(patch_id_map,f))==ndaf.end() ) 
         {
            ndaf[get(patch_id_map,f)] = {tri};
         }
         else 
         {
            ndaf[get(patch_id_map,f)].push_back(tri);
         }
         
     }
     
     // sort on second element 
     return ndaf;
     //return Tri2tagvec;
   }
      
  // DocString: remove_small_components   
  /** 
   * @brief Removes connected components whose area or volume is under a certain threshold value. 
   * @precondition surface mesh must bound a volume.
   * @param volume_threshold ratio value of the volume such that only connected components whose volume is larger than this value are kept (only applies to closed connected components) 
   * @see[remove_connected_components_of_negligible_size](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__repairing__grp.html#gac544fcaba1d59d330a3a1536caff392a)
   */   
   void remove_small_components(double volume_threshold)
   {
       if( !does_bound_a_volume() ) 
            throw PreconditionError("Surface does not bound a volume.");
       CGAL::Polygon_mesh_processing::remove_connected_components_of_negligible_size(mesh, CGAL::parameters::volume_threshold( volume_threshold*volume()));
   }

  /**
   * @brief checks vertices position property relative to surface.
   *
   * This function is used to investigate if a surface is partially 
   * inside and outside of another surface.
   *  
   * @tparam CGAL::Bounded_side. 
   * @see [CGAL::Bounded_side] (https://doc.cgal.org/latest/Kernel_23/group__kernel__enums.html)
   *
   * @param other SVMTK Surface object. 
   * @returns  first true if property match first template, second true if property match second template.
   */ 
   template< CGAL::Bounded_side A , CGAL::Bounded_side B>
   std::pair<bool,bool> check_vertices(Surface &other)
   {
      assert_non_empty_mesh();
      bool query1 = false, query2 = false; 
      Inside is_inside_query(other.get_mesh()); 
      
      
      for(vertex_descriptor vit : mesh.vertices() )
      {
         CGAL::Bounded_side res =  is_inside_query(mesh.point(vit));
         if( res==A )                                
            query1 = true;
         if( res==B )
            query2 = true;
      }
      return std::make_pair(query1,query2);
   }

   template< CGAL::Bounded_side A , CGAL::Bounded_side B>
   std::pair<bool,bool> check_points( point_vector points)
   {
      assert_non_empty_mesh();
      bool query1 = false, query2 = false; 
      Inside is_inside_query( this->mesh); 
      
      for(auto pit : points )
      {
         CGAL::Bounded_side res =  is_inside_query(pit);
         if( res==A )                                
            query1 = true;
         if( res==B )
            query2 = true;
      }
      return std::make_pair(query1,query2);
   }
   
 
   // DocString: repair_self_intersections  
  /**
   * @brief Removes self intersection from surface.
   * @note This function uses experimental CGAL algorithms to remove self-intersections.
   * @see [remove_self_intersections](https://github.com/CGAL/cgal/blob/master/Polygon_mesh_processing/include/CGAL/Polygon_mesh_processing/repair_self_intersections.h)
   *     
   * @param volume_threshold ratio value of the volume such that only connected components whose volume is larger than this value are kept 
            (only applies to closed connected components) 
   * @param cap_thresold anlge between 160 180 degrees 
   *         CGAL explanation->the cosine of a minimum angle such that if `f` has an angle greater than this bound,
   *                             it is a cap. The threshold is in range `[-1 0]` and corresponds to an angle 
   *                             between `90` and `180` degrees.
   * @param needle_threshold long edge divided by short edge, i.e ratio  
   * @param collapse_threshold ratio of the average edge length. 
   * @returns std::pair<bool,int>  true if complete, and number of remaining self intersections.
   * # TODO Investigate crash 
   */
   std::pair<bool,int> repair_self_intersections(double volume_threshold=0.01, double cap_threshold=170, double needle_threshold=2.7, double collapse_threshold=0.14)
   {
        double avgel = average_edge_length();
        CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
        CGAL::Polygon_mesh_processing::remove_connected_components_of_negligible_size(mesh,CGAL::parameters::volume_threshold( volume_threshold*volume()));
        CGAL::Polygon_mesh_processing::remove_degenerate_edges(edges(mesh), mesh ); 
        CGAL::Polygon_mesh_processing::remove_almost_degenerate_faces( faces(mesh), mesh, CGAL::parameters::needle_threshold( needle_threshold).
                                               cap_threshold(cap_threshold).
                                               collapse_length_threshold(avgel*collapse_threshold));
                                             
        CGAL::Polygon_mesh_processing::remove_connected_components_of_negligible_size(mesh,CGAL::parameters::volume_threshold( volume_threshold*volume()));     
        CGAL::Polygon_mesh_processing::experimental::remove_self_intersections(mesh);

        auto intersecting_faces = get_self_intersecting_faces();
        this->remove_faces(intersecting_faces);
        int after = num_self_intersections();
        return std::make_pair(after==0, after); 
   }

  /**
   * @brief Returns the vertices of the surface mesh
   * @returns a vector of surface vertices 
   */
   vertex_vector get_vertices() 
   {
      assert_non_empty_mesh();
      vertex_vector result;
      for(vertex_descriptor vit : mesh.vertices())
          result.push_back(vit);
      return result;   
   }

  /**
   * @brief Finds and returns vertices that are inside another SVMTK Surface object.
   * @tparam CGAL::Bounded_side A.
   * @tparam CGAL::Bounded_side B.
   * @see [CGAL::Bounded_side] (https://doc.cgal.org/latest/Kernel_23/group__kernel__enums.html)
   * @param[in] other SVMTK Surface object. 
   * @param[in,out] vertices vector of this->surface mesh vertices.   
   * @returns vector of vertices. 
   * @overload 
   */
   template< CGAL::Bounded_side A , CGAL::Bounded_side B>
   vertex_vector_map get_vertices_with_property(Surface &other,vertex_vector_map &mvertices)
   {
      assert_non_empty_mesh();
      Inside is_inside_query(other.get_mesh()); 
      for(auto vit = mvertices.begin(); vit!= mvertices.end(); )
      {
         CGAL::Bounded_side res = is_inside_query(mesh.point(vit->first));
         if( res==A or res==B )
           vit++;
         else 
           mvertices.erase(vit++);       
      }
      return mvertices;
   }
   
   /*bool remove_close_vertices(vertex_vector &mvertices)
   {
         FT avg_edge_length = FT(average_edge_length());  
         vertex_vector to_remove;
         for( auto vit = mvertices.begin(); vit!= mvertices.end(); vit++)
         {
              CGAL::Halfedge_around_target_circulator<Mesh> vbegin(mesh.halfedge(*hit), mesh), done(hbegin);
              do
              {
                 edge   = Vector_3(current,mesh.point(*vbegin));
                 edge_length = CGAL::sqrt(edge.squared_length());  
                 edge_length = CGAL::min(avg_edge_length, edge_length);
                 if( edge_length <= 0.1*avg_edge_length ) 
    
                *hbegin++;
              }while(hbegin!=done);
         }
         
         for( auto vit : to_remove) 
         {
            CGAL::
         
         
         }
   
   
   }*/
   
   
   
  /**
   * @brief Finds and returns vertices that are inside another SVMTK Surface object.
   * @tparam CGAL::Bounded_side A.
   * @tparam CGAL::Bounded_side B.
   * @see [CGAL::Bounded_side] (https://doc.cgal.org/latest/Kernel_23/group__kernel__enums.html)
   * @param other SVMTK Surface object. 
   * @param vertices vector of this-> surface mesh vertices.  
   * @returns a vector of vertices. 
   * @overload 
   */
   template< CGAL::Bounded_side A , CGAL::Bounded_side B>
   vertex_vector get_vertices_with_property(Surface &other,vertex_vector &vertices)
   {
      assert_non_empty_mesh();
      Inside is_inside_query(other.get_mesh()); 
      for(auto vit = vertices.begin(); vit!= vertices.end(); )
      {
         CGAL::Bounded_side res =  is_inside_query(mesh.point(*vit));
         if( res==A or res==B )
           vit++;
         else 
           vertices.erase(vit);       
      }
      return vertices;
   }
   
  /**
   * @brief Finds and returns vertices that are inside another SVMTK Surface object.
   * @tparam CGAL::Bounded_side A.
   * @tparam CGAL::Bounded_side B.
   * @see [CGAL::Bounded_side] (https://doc.cgal.org/latest/Kernel_23/group__kernel__enums.html)
   * @param other SVMTK Surface object.
   * @returns a vector of vertices. 
   * @overload 
   */
   template<CGAL::Bounded_side A, CGAL::Bounded_side B>
   vertex_vector get_vertices_with_property(Surface &other)
   {
     vertex_vector vertices = this->get_vertices();
     return get_vertices_with_property<A,B>(other,vertices);
   }   

  /** 
   * @brief Finds and returns vertices that are outside another SVMTK Surface object..
   * @param other SVMTK Surface object.
   * @param vertices vector of surface mesh vertices.    
   * @returns a vector of vertices. 
   * @overload 
   */
   template<typename Vertex_iterator>
   Vertex_iterator get_vertices_outside(Surface &other,Vertex_iterator &vertices) 
   {
      return this->get_vertices_with_property<CGAL::ON_UNBOUNDED_SIDE,CGAL::ON_UNBOUNDED_SIDE>(other, vertices);
   }
   
   /**  
    * @brief Finds and returns vertices that are outside another SVMTK Surface object.
    * @param other SVMTK Surface object.
    * @returns a vector of vertices. 
    */
    vertex_vector get_vertices_outside(Surface &other) 
    {
        vertex_vector vertices = this->get_vertices();
        return this->get_vertices_with_property<CGAL::ON_UNBOUNDED_SIDE,CGAL::ON_UNBOUNDED_SIDE>(other, vertices);
    }

   /** 
    * @brief Finds and returns vertices that are inside another SVMTK Surface object.
    * @param other SVMTK Surface object.
    * @param vertices vector of surface mesh vertices.   
    * @returns a vector of vertices.
    * @overload     
    */
    template<typename Vertex_iterator> 
    Vertex_iterator get_vertices_inside(Surface &other,Vertex_iterator &vertices)  
    {
       return this->get_vertices_with_property<CGAL::ON_BOUNDED_SIDE,CGAL::ON_BOUNDARY>(other, vertices);
    }

  /** 
   * @brief Finds and returns vertices that are on the boundary of another SVMTK Surface object.
   * @param other SVMTK Surface object.
   * @param vertices vector of surface mesh vertices.   
   * @returns a vector of vertices. 
   */
   template<typename Vertex_iterator>  
   Vertex_iterator get_vertices_on_boundary(Surface &other,Vertex_iterator &vertices) 
   {
       return this->get_vertices_with_property<CGAL::ON_BOUNDARY,CGAL::ON_BOUNDARY>(other, vertices);
   }

   /** 
    * @brief Finds and returns vertices that are inside another SVMTK Surface object.
    * @param other SVMTK Surface object.
    * @returns a vector of vertices.
    * @overload 
    */
    vertex_vector get_vertices_inside(Surface &other) 
    {
       vertex_vector vertices = this->get_vertices();
       return this->get_vertices_with_property<CGAL::ON_BOUNDED_SIDE,CGAL::ON_BOUNDARY>(other,vertices);
    }

   /**
    * @brief Finds and returns surface mesh vertices that are close to another SVMTK Surface object.
    * @tparam mode an integer that decides 
    * @param other SVMTK Surface object.
    * @returns a vector of vertices. 
    */
    vertex_vector get_close_vertices(Surface &other)
    {
        return this->get_close_vertices_with_property<0,false>(other);
    }    
  
   /**
    * @brief Finds and returns surface mesh vertices that are close to another SVMTK Surface object.
    * @param other SVMTK Surface object.
    * @returns a vector of vertices. 
    */
    vertex_vector get_close_vertices(Surface &other, vertex_vector &mvertices)
    {
       return this->get_close_vertices_with_property<0,false>(other,mvertices);
    } 


   /**
    * @brief Finds close vertices where the direction can be specificed as either 
    *        negative normal direction or positive normal direction. 
    *
    * @tparam A  sets the first element of the  search. If other is (*this), 
    *            then A=1 so that search vertex is not chosen as the closest vertex.  
    * @tparam B If true; close vertices in negative normal direction are skipped.
                If false; close vertices in positiv normal direction are skipped.
    * @param other SVMTK surface object, itself is an option.
    * @returns a vector of vertices. 
    */
    template<int A, bool B> 
    vertex_vector get_close_vertices_with_property(Surface &other)
    {
       vertex_vector vertices =this->get_vertices();
       return get_close_vertices_with_property<A,B>(other, vertices);
    }  
                    
   /** 
    * @brief Separates non-connected surface mesh vertices.
    * Separates non-connected surface mesh vertices so that the nearest surface mesh vertices are adjacent.  
    * @param adjustment multiplier of the edge movement.
    * @returns the number of adjusted vertices.  
    */
    vertex_vector get_narrow_gaps() 
    {
       vertex_vector vertices = this->get_vertices();
       return this->get_close_vertices_with_property<1,true>(*this,vertices);
    }
      
   /**
    * @brief Finds and returns non-adjacent vertices that are close in a negative normal direction.
    * @param  other SVMTK Surface object
    * @returns a vector of vertices that are close in a negative normal direction.
    */
    vertex_vector get_needle_vertices()
    {
       vertex_vector vertices = this->get_vertices();
       return this->get_close_vertices_with_property<1,false>(*this, vertices);
    }
   

     
   /**
    * @brief Finds close vertices where the direction can be specificed as either 
    *        negative normal direction or positive normal direction. 
    *
    * @tparam A  sets the first element of the  search. If other is (*this), 
    *            then A=1 so that search vertex is not chosen as the closest vertex.  
    * @tparam B If true; close vertices in negative normal direction are skipped.
                If false; close vertices in positiv normal direction are skipped.
    * @param other SVMTK surface object, itself is an option.
    * @param[in,out] mvertices a vector of vertices. 
    * @returns a vector of vertices. 

    */
    template<int A, bool B> 
    vertex_vector get_close_vertices_with_property(Surface &other, vertex_vector& mvertices) 
    {                     
       assert_non_empty_mesh();
       CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);

       vertex_vector results;
       Vertex_point_pmap vppmap = get(CGAL::vertex_point,other.get_mesh());  

       Tree tree(vertices(other.get_mesh()).begin(),
                 vertices(other.get_mesh()).end(),
                 Splitter(),
                 Traits(vppmap));
 
       Distance tr_dist(vppmap);
       FT distance, edge_length;
       Point_3 closest;
       Vector_3 direction,normal,edge;
       bool flag;

       FT avg_edge_length = FT(average_edge_length());  

       for(auto vit = mvertices.begin(); vit!= mvertices.end(); )
       {
          K_neighbor_search search(tree, mesh.point(*vit), 2,0,true,tr_dist); 
          closest = other.get_mesh().point((search.begin()+A)->first);
          Point_3 current = mesh.point(*vit);
          direction =  Vector_3(current,closest);
          
          normal =  CGAL::Polygon_mesh_processing::compute_vertex_normal(*vit,mesh);
          distance = CGAL::sqrt(direction.squared_length());

          // Just for separate narrow gaps, so that needles are excluded
          if( (normal*direction<=0)==true and B) 
          { 
            mvertices.erase(vit);
            continue;
          }
          flag = true;
          CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(*vit), mesh), done(vbegin);
          do
          {
             edge   = Vector_3(current,mesh.point(*vbegin));
             edge_length = CGAL::sqrt(edge.squared_length());  
             
             edge_length = CGAL::min(avg_edge_length, edge_length);
             if( _prratio*edge_length <= distance or closest == mesh.point(*vbegin) ) 
             {
               flag=false;
               break;
             }
             *vbegin++;
          }while(vbegin!=done);

          if (flag)
             vit++;
          else 
             mvertices.erase(vit);
       }
       return mvertices;
    }  
       

  // DocString: separate_narrow_gaps
   /**
    * @brief Separates unconnected surface mesh vertices that are close in a positive normal direction by
    *        moving vertices in a negative surface normal direction, i.e. contraction.
    *  
    * @param adjustment multiplier of the edge movement.
    * @returns std::pair<bool,int> completion before maximum iteration reached, and number of adjusted vertices.  
    */  
    std::pair<bool,int> separate_narrow_gaps(double adjustment, double smoothing, int max_iter)  
    {
       assert_non_empty_mesh();
       CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh); 
       return manipulate_close_vertex_selection<1,true>(*this, -abs(adjustment),smoothing,max_iter);
    }

  // DocString: separate_close_vertices
  /** 
   * @brief Separates close non-connected surface mesh vertices.
   * 
   * Separates close non-connected surface mesh vertices so that the closest surface mesh vertex is one
   * that is connected, i.e. shares a edge.
   * 
   * @param adjustment multiplier for vertex movement.
   * @param max_iter maximum number of iterations for a while loop.
   * @returns std::pair<bool,int> true if completed and number of adjusted vertices. 
   */
   std::pair<bool,int> separate_close_vertices(double adjustment, int max_iter) 
   {
      assert_non_empty_mesh();
      double smoothing = adjustment;
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh); 
      if( adjustment<0 )
         throw  InvalidArgumentError("Adjusment must be positive.");
      return manipulate_close_vertex_selection_with_direction<1,false>(*this,adjustment, smoothing ,max_iter); 
   }

  // DocString: embed
  /**
   * @brief Moves all vertices so that they are inside another surface. 
   *
   * @param other SVMTK Surface object. 
   * @param adjustment multiplier for vertex movement.
   * @param smoothing Laplacian smoothing factor for moved vertices after each iteration in while loop. @see laplacian_smoothing
   * @param max_iter maximum number of iterations for a while loop.
   *
   * @returns std::pair<bool,int> true if completed and number of adjusted vertices. 
   * @throws InvalidArgumentError if adjustment is positive, i.e. expansion instead of contraction.
   */
   std::pair<bool,int> embed(Surface& other, double adjustment,  double smoothing, int max_iter) 
   {
      assert_non_empty_mesh();
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh); 
      return manipulate_vertex_selection<CGAL::ON_UNBOUNDED_SIDE,CGAL::ON_BOUNDARY>(other, -abs(adjustment), smoothing, max_iter);  
   }
 
  // DocString: expose
  /**
   * @brief Moves vertices in a negative normal direction so that they are outside specified surface.
   *
   * @param other SVMTK Surface object. 
   * @param adjustment multiplier for vertex movement.
   * @param smoothing Laplacian smoothing factor for moved vertices after each iteration in while loop. @see laplacian_smoothing
   * @param max_iter maximum number of iterations for a while loop.
   * @returns std::pair<bool,int> true if completed and number of adjusted vertices. 
   * @throws InvalidArgumentError if adjustment is positive, i.e. expansion instead of contraction.
   */
   std::pair<bool,int> expose(Surface& other, double adjustment,  double smoothing, int max_iter) 
   {
      assert_non_empty_mesh();
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh); 
     return manipulate_vertex_selection<CGAL::ON_BOUNDED_SIDE,CGAL::ON_BOUNDARY>(other, -abs(adjustment), smoothing, max_iter);
   }
  
  // DocString: enclose
  /** 
   * @brief Moves all vertices so that they are outside another surface.  
   * 
   Precondition that the surfaces intersect. 
   * @param other SVMTK Surface object.
   * @param adjustment multiplier for vertex movement.
   * @param smoothing Laplacian smoothing factor for moved vertices after each iteration in while loop. @see laplacian_smoothing
   * @param max_iter maximum number of iterations for a while loop.
   * @returns std::pair<bool,int> true if completed and number of adjusted vertices. 
   * @throws InvalidArgumentError if adjustment is negative, i.e. contraction instead of expansion.
   */
   std::pair<bool,int> enclose(Surface& other, double adjustment ,double smoothing , int max_iter) 
   {
      assert_non_empty_mesh();
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh); 
      return manipulate_vertex_selection<CGAL::ON_BOUNDED_SIDE,CGAL::ON_BOUNDARY>(other, abs(adjustment), smoothing, max_iter);   
   }

   // DocString: separate 
  /**
   * @brief  Moves all vertices that are close to another surface in opposite direction.
   * @precondition the other surface must either enclosed or be embedded in the Surface object. 
   * @param other SVMTK Surface object.
   * @param adjustment multiplier for vertex movement.
   * @param max_iter maximum number of iterations for a while loop.
   * @returns std::pair<bool,int> true if completed and number of adjusted vertices. 
   * @throws InvalidArgumentError if surfaces intersect.
   */
   std::pair<bool,int> separate(Surface& other, double adjustment, double smoothing, int max_iter)
   { 
      assert_non_empty_mesh();
      
      //std::pair<bool,bool> queries = other.check_vertices<CGAL::ON_BOUNDED_SIDE,CGAL::ON_UNBOUNDED_SIDE>(*this);
      //queries = check_vertices<CGAL::ON_BOUNDED_SIDE,CGAL::ON_UNBOUNDED_SIDE>(other);
      return manipulate_close_vertex_selection<0,false>(other, adjustment, smoothing, max_iter);

   }

   /**
    * @brief Moves the vertex position of a selection of Surface vertices dependent on different queries, like inside or outside another surface.
    * @tparam CGAL::Bounded_side.
    * @param other SVMTK Surface object.
    * @param adjustment times the shortest connected edge length gives the limit to the vertex movement. 
    * @param smoothing laplacian smoothing parameter for smoothing after each iteration. 
    * @param max_iter the maximum number of iteration for the while. 
    * @returns std::pair<bool,int>. True if algorithm complets before max_iter is reached. The integer is the number of vertex manipulation that
    * is completed.  
    */
    template< CGAL::Bounded_side A , CGAL::Bounded_side B> 
    std::pair<bool,int> manipulate_vertex_selection(Surface &other, double adjustment,  double smoothing, int max_iter)
    {

       double smth = smoothing;
       vertex_vector vertices = get_vertices_with_property<A,B>(other);
        
       double displacement_bound = get_lower_displacment_bound();
       double srf = get_smoothing_reduction_factor();

       int iter = 0;
       while( !vertices.empty() )
       {     
          // Relative to the surface normal, i.e. negative adjusment gives implosion 
          vertex_vector_map vertex_displacement = get_vertex_displacement(vertices, adjustment, smth); 

          auto max_element = std::max_element( vertex_displacement.begin(), vertex_displacement.end(), find_largest_displacement());    
          
          set_adjacent_vertices(vertex_displacement);
          
          if ( max_element->second.squared_length() < displacement_bound ) 
               smth = srf*smth;

          vertices = get_vertices_with_property<A,B>(other,vertices);
          
          if( this->num_self_intersections()>0 ) 
          {
             this->repair_self_intersections(); 
             vertices = get_vertices_with_property<A,B>(other);
             smth = smoothing;     
          }
          
           if( vertices.empty() ) 
           {
              vertices = get_vertices_with_property<A,B>(other);
              smth = smoothing;
           }
                       
          if( ++iter>max_iter )
             return std::make_pair(false, vertices.size() );
       }
       // 
       //manipulate_close_vertex_selection<0,false>(other,adjustment,smoothing,10);
       //vertices = get_vertices_with_property<A,B>(other,vertices);
       return std::make_pair(true, vertices.size() );
    }
      
   /**
    * @brief Moves the vertex position of a selection of Surface vertices dependent 
    * on adjaceny to other vertices, itself and other surfaces.
    *
    * @tparam int indicating the starting point of the adjacency search. If other is (*this) then A=1 so that vertex does not choose itself as a 
    *             close vertex.  
    * 
    * @param other SVMTK Surface object.
    * @param adjustment times the shortest connected edge length gives the limit to the vertex movement. 
    * @param smoothing laplacian smoothing parameter for smoothing after each iteration. 
    * @param max_iter the maximum number of iteration for the while. 
    *
    * @returns std::pair<bool,int>. True if algorithm complets before max_iter is 
    * reached. The integer is the number of vertex manipulation that
    * is completed.  
    */
    template<int A, bool B = false>
    std::pair<bool,int> manipulate_close_vertex_selection(Surface &other, double adjustment,  double smoothing, int max_iter)
    { 
       double displacement_bound = get_lower_displacment_bound();
       double srf = get_smoothing_reduction_factor();

       double smth = smoothing;
       vertex_vector vertices = get_close_vertices_with_property<A,B>(other);

       int iter = 0;


       while( !vertices.empty() )
       {       
          vertex_vector_map vertex_displacement = get_vertex_displacement(vertices, adjustment, smth);
          
          auto max_element = std::max_element( vertex_displacement.begin(), vertex_displacement.end(), find_largest_displacement());    
         
          set_adjacent_vertices(vertex_displacement);
          
          if ( max_element->second.squared_length() < displacement_bound) 
               smth = srf*smth;
                   
          vertices = get_close_vertices_with_property<A,B>(other,vertices); 

          if( this->num_self_intersections()>0 ) 
          {
             this->repair_self_intersections(); 
             vertices = get_close_vertices_with_property<A,B>(other);
             smth= smoothing;    
          }
          if( vertices.empty()  )  
          {
              vertices = get_close_vertices_with_property<A,B>(other); 
              smth = smoothing;
          }          
          if( ++iter>max_iter )
             return std::make_pair(false, vertices.size() );
       }
       return std::make_pair(true,  vertices.size() );
    }

   /**  
    * Experimental:
    * @brief Moves the vertex position of Surface vertices that are close to none-connected  
    *        vertices.
    *
    * @tparam A sets the first element of the adjacency search. If other is (*this), 
    *           then A=1 so that search vertex is not chosen as the closest vertex.  
    * 
    * @param other SVMTK Surface object.
    * @param adjustment multipler of the shortest connected edge length, which gives the length of the vertex displacement.
    * @param smoothing laplacian smoothing parameter for smoothing after each iteration. 
    * @param max_iter the maximum number of iteration for the while. 
    *
    * @returns std::pair<bool,int>. True if algorithm complets before max_iter is reached. The integer is
    *         the number of vertex manipulation that is completed.  
    *  
    */
    template<int A, bool B>
    std::pair<bool,int> manipulate_close_vertex_selection_with_direction(Surface &other, double adjustment, double smoothing,int max_iter)
    {   
       double displacement_bound = get_lower_displacment_bound();
       double srf = get_smoothing_reduction_factor();

       double smth = smoothing;  
       vertex_vector vertices  =  get_close_vertices_with_property<A,B>(other); 
       add_adjacent_vertices(vertices,60);


       int iter = 0;

       while( !vertices.empty() )
       {  
          // adjustment is in negative direction 
          vertex_vector_map vertex_displacement = get_close_vertex_displacement<A>(other, vertices, adjustment, smoothing); 
          
          auto max_element = std::max_element( vertex_displacement.begin(), vertex_displacement.end(), find_largest_displacement());    
          
          if ( max_element->second.squared_length() < displacement_bound) 
               smth = srf*smth;
        
          set_adjacent_vertices(vertex_displacement);
          
          vertices = get_close_vertices_with_property<A,B>(other, vertices); 

          if( this->num_self_intersections()>0 ) 
          {
             this->repair_self_intersections(); 
             vertices = get_close_vertices_with_property<A,B>(other);
             add_adjacent_vertices(vertices,60);                
             smth = smoothing;   
          }
               
          if( vertices.empty()  )  
          {   
              vertices  = get_close_vertices_with_property<A,B>(other);
              add_adjacent_vertices(vertices,60);       
              smth = smoothing;   
          }
                    
          if( ++iter>max_iter )
             return std::make_pair(false, vertices.size() );
       }
             
       return std::make_pair(true,  vertices.size() );
    }
 

  /**
   * @brief Checks if the surface bounds a volume.
   * @returns true if surface bounds a volume  
   */
   bool does_bound_a_volume()
   {  
      return CGAL::is_closed(mesh);
   }

  // DocString: distance_to_point
  /**
   * @brief Computes the distance between a point and the closest vertex in the triangulated surface.
   * @param point SVMTK Surface Point_3 object. 
   * @returns the distance from point and the closest vertex on the surface
   */
   double distance_to_point(Point_3 point)
   {
      assert_non_empty_mesh();
      vertex_vector vertices = get_closest_vertices(point,1);
      if( vertices.size()==0 )
         return 0;
      Vector_3 vector(mesh.point(vertices[0]),point); 
      return CGAL::to_double(CGAL::sqrt(vector.squared_length()));
   }

  /**
   * @brief Handles the output for PreconditionError in boolean operations.
   * @param SVMTK Surface object.
   */
   std::string CGAL_precondition_evaluation(Surface other)
   {
      std::string output = "Following preconditions failed: ";       
      if (CGAL::Polygon_mesh_processing::does_bound_a_volume(mesh))
         output = output +"\n" +"- Surface does no bound a volume.";       
      if (CGAL::Polygon_mesh_processing::does_bound_a_volume(other.get_mesh()))
         output = output +"\n" + "- Argument surface does no bound a volume.";  
      if ( CGAL::Polygon_mesh_processing::does_self_intersect(mesh)) 
         output = output +"\n" + "- Surface self intersection in witihn union volume.";      
      if ( CGAL::Polygon_mesh_processing::does_self_intersect(other.get_mesh())) 
         output = output +"\n" + "- Argument surface self intersection in witihn union volume.";
      return output;
   }

  /**
   * @brief Sets the face oritentiaton outwards 
   */
   void set_outward_face_orientation()
   {
     assert_non_empty_mesh();
     if( CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)) )
        CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
   }

  // DocString: intersection
  /**
   * @brief Computes the intersection between two triangulated surface mesh.
   * @see [corefine_and_compute_intersection](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html) 
   * The functions have precondition that both surfaces bounds a volume  and that both surfaces does not have self-intersections
   * where the surfaces intersect.
   * @param SVMTK Surface object
   * @returns true if intersection computation is successful  
   */
   bool surface_intersection(Surface other)
   {
      assert_non_empty_mesh();
      other.assert_non_empty_mesh();
      try
      {
        return CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh, other.get_mesh(),mesh);
      }
      catch(const std::exception &exc)
      {
        CGAL_precondition_evaluation( other);
        std::string output = "CGAL precondition error\n"+ CGAL_precondition_evaluation(other);
        throw PreconditionError(output.c_str());
      }       
   }  
   
  // DocString: difference
  /**
   * @brief Computes the difference between two triangulated surface mesh.
   * @see [corefine_and_compute_difference](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html) 
   * The functions have precondition that both surfaces bounds a volume and that both surfaces does not have self-intersections
   * where the surfaces intersect.
   * @param other SVMTK Surface object
   * @returns true if intersection computation is successful  
   */
   bool surface_difference(Surface other)
   {    
     assert_non_empty_mesh();
     other.assert_non_empty_mesh();
     try 
     {
       return CGAL::Polygon_mesh_processing::corefine_and_compute_difference(mesh , other.get_mesh(), mesh);
     }
     catch(const std::exception &exc)
     {
        CGAL_precondition_evaluation(other);
        std::string output = "CGAL precondition error\n"+ CGAL_precondition_evaluation(other);
        throw PreconditionError(output.c_str());
     }   
   }

  // DocString: union
  /**
   * @brief Computes the union between two triangulated surface mesh.
   * @see [corefine_and_compute_union](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html) 
   * The functions have precondition that both surfaces bounds a volume and that both surfaces does not have self-intersections 
   * where the surfaces intersect.
   * @param other SVMTK Surface object
   * @returns true if intersection computation is successful  
   */
   bool surface_union(Surface other)
   {  
     assert_non_empty_mesh();
     other.assert_non_empty_mesh();
     if( !does_bound_a_volume()) 
         return false;
     if( !other.does_bound_a_volume())
         return false;
     try 
     { 
       return CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh, other.get_mesh(), mesh);
     }
     catch (const std::exception &exc)
     {
       std::string output = "CGAL precondition error\n"+ CGAL_precondition_evaluation(other); 
       throw PreconditionError(output.c_str());
     }   
   }

  // DocString: keep_largest_connected_component  
  /**
   * @brief Keeps the largest connected component of the stored mesh.
   * @see [keep_largest_connected_components](https://doc.cgal.org/latest/Polygon_mesh_processing/group__keep__connected__components__grp.html) 
   * @returns the number of connected components.
   */ 
   int keep_largest_connected_component()
   {      
       assert_non_empty_mesh();
       CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
       return CGAL::Polygon_mesh_processing::keep_largest_connected_components(mesh,1);
   }
  
  /** 
   * 
   * 
   */
   int num_connected_component()
   {      
       assert_non_empty_mesh();
       CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
       // Simplify 
       Mesh::Property_map<face_descriptor, std::size_t> fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
       return CGAL::Polygon_mesh_processing::connected_components(mesh,fccmap);
   }

   
   std::vector< std::shared_ptr<Surface> > connected_components() 
   {
       assert_non_empty_mesh();
       CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
       std::vector<Mesh> cc_meshes; 
       CGAL::Polygon_mesh_processing::split_connected_components(mesh, cc_meshes );
       std::vector< std::shared_ptr<Surface> > result; 
       for( auto cc_mesh : cc_meshes) 
       {
           result.emplace_back(std::shared_ptr<Surface>(new Surface(cc_mesh)));
       }
       return result;
   }
   




  /**
   * @brief Finds a specified number mesh vertices that are closest to a point outside the mesh. 
   * @param p1 is an arbiraty point outside the surface mesh 
   * @param num the maximum number of vertices that are added to the result
   * @returns a vertex_vector that holds sum of `values`, or 0.0 if `values` is empty.
   */
   vertex_vector get_closest_vertices(Point_3 p1, int num=8)
   {
     Vertex_point_pmap vppmap = get(CGAL::vertex_point,mesh); 
     vertex_vector results;

     Tree tree(vertices(mesh).begin(),
               vertices(mesh).end(),
               Splitter(),
               Traits(vppmap));

     Distance tr_dist(vppmap);
     K_neighbor_search search(tree, p1, num, 0, true, tr_dist); 
 
     for(auto vit=search.begin(); vit!=search.end(); ++vit) 
        results.push_back(vit->first); 
     return results;
   } 

  // DocString: get_closest_points
  /**
   * @brief Finds a specified number mesh points that are closest to a point not on the mesh. 
   * @param point is an arbiraty point not on the surface mesh 
   * @param num the maximum number of points that are added to the result
   * @returns a point_vector that holds sum of `values`, or 0.0 if `values` is empty.
   */
   point_vector get_closest_points(Point_3 point, int num=8)
   {
      vertex_vector vertices = get_closest_vertices( point, num);
      return get_points(vertices);
   } 




  // DocString: cylindrical_extension
  /** 
   * @brief Constructs a cylindrical extension, which is a cylinder surface mesh combined with a sphere 
   * on the end closest to the mesh.
   * 
   * The user can use the union operation to combine it with the main mesh. 
   * The cylinder surface mesh is determined by centeriod of vertex points that are closest to a point 
   * outisde the mesh, radius, length and the option of crating a cylinder normal to the surface mesh. 
   * Works best on convex surfaces. 
   *
   * @param point is an arbiraty point outside the surface mesh, the closest vertex is 
   * @param radius of the cylinder surface
   * @param length of the cylinder surface from the surface mesh.
   * @param edge_length the edge length of the output surface.
   * @param normal if true the cylinder is set normal to the surface mesh, else in the direction of point. 
   * @returns a SVMTK Surface object.
   * @throws InvalidArgumentError if point is inside surface.
   */
   std::shared_ptr<Surface> cylindrical_extension(const Point_3& point, double radius, double length, double edge_length, bool use_normal=true)
   {
      assert_non_empty_mesh();
 
      if( length<=0 )
         throw InvalidArgumentError("The length must be larger than 0.");         
           
      Vector_3 normal;
     
      vertex_vector vertices = get_closest_vertices(point,6);
      
      std::vector<Point_3> points = get_points(vertices);  
      
      Point_3 closest_point;
         
      if (use_normal) 
      {
          for(auto i : vertices ) 
          {
              normal += CGAL::Polygon_mesh_processing::compute_vertex_normal(i,mesh);
          }
          if ( CGAL::sqrt(normal.squared_length()) >0)
          {
             normal /= CGAL::sqrt(normal.squared_length());
          }
          else 
          {
             normal = CGAL::Polygon_mesh_processing::compute_vertex_normal(vertices[0],mesh);
          }
          closest_point = CGAL::centroid(points.begin(),points.end(),CGAL::Dimension_tag<0>()); 
      }
      else 
      {
          closest_point = mesh.point(vertices[0]); 
                    
          if( point==closest_point)   
          {
             normal = CGAL::Polygon_mesh_processing::compute_vertex_normal(vertices[0],mesh);
          }
          else 
          {
             normal = Vector_3(closest_point, point);
          }
      }
      
      Point_3 target = closest_point + normal*length;  
      Point_3 source = closest_point - normal*radius;
      
      std::shared_ptr<Surface> cylinder(new Surface()); 
      
      cylinder->make_cylinder(source,target,radius,edge_length); 
      return cylinder; 
   }
   
  // DocString: cylindrical_extension
  /**
   * @brief Creates a triangulates cylinder surface from a point outside the surface to the closest surface point with given radius, 
   *  length of cylinder and edge length. The option of normal can be used for the cylinder to be perpendicular to the surface.     
   * 
   * @param x coordinate of point.
   * @param y coordinate of point.  
   * @param z coordinate of point. 
   * @param radius of the cylinder.
   * @param length of the cylinder.
   * @param edge length of the triangles.
   * @param normal if true the cylinder centerline will be normal on the surface.
   * @returns a SVMTK Surface object.
   * @overload 
   */   
   std::shared_ptr<Surface> cylindrical_extension(double x, double y ,double z, double radius, double length,double edge_length, bool normal=true)                
   {
     return cylindrical_extension(Point_3(x,y,z),radius,length,edge_length,normal);
   }
   
  // DocString: cylindrical_connection
  /**
   * @brief Constructs a cylindrical connection bridge between the shortests line between two points 
   * in different surface.   
   *
   * @note The user can use the union operation to combine it with the main mesh. 
   * @param other SVMTK Surface object.
   * @param radius of the cylinder surface
   * @param edge_length the edge length of the output surface.
   * @returns a SVMTK Surface object.
   */
   std::shared_ptr<Surface> cylindrical_connection(Surface other, double radius, double edge_length)
   {
      auto points = find_closest_point_pair(other);
      
      Vector_3 direction = Vector_3(points.first , points.second);
       
      FT length = CGAL::sqrt(direction.squared_length());
      direction/=length;  
 
      Point_3 p2 = points.first  - direction*radius; 
      Point_3 p3 = points.second + direction*radius;
   
      std::shared_ptr<Surface> cylinder(new Surface()); 
      cylinder->make_cylinder(p2,p3,radius,edge_length); 
      return cylinder;
   }
   
  /**
   * @brief Finds the closest pair of points between *this and another 
   *  SVMTK Surface object. 
   * 
   * @param other SVMTK Surface object.
   * @return a pair of SVMTK Point_3 objects  
   */
   std::pair<Point_3,Point_3> find_closest_point_pair(Surface other) 
   {
      Vertex_point_pmap vppmap = get(CGAL::vertex_point,other.get_mesh());  

      Tree tree(vertices(other.get_mesh()).begin(),
                vertices(other.get_mesh()).end(),
                Splitter(),
               Traits(vppmap)
      );
      Distance tr_dist(vppmap);
      FT distance;   
      
      Point_3 p0 = this->get_points()[0];
      Point_3 p1 = other.get_points()[0];
      
      FT min_distance = FT( CGAL::squared_distance(p0,p1) ); // std::numeric_limits<double>::max() 
      Point_3 query_point;

      Vector_3 direction;
      
      for(vertex_descriptor vit : mesh.vertices())
      {
         K_neighbor_search search(tree, mesh.point(vit), 2,0,true,tr_dist); 
         
         query_point = other.get_mesh().point((search.begin())->first);
         Point_3 current = mesh.point(vit);
         
         distance = CGAL::squared_distance(current,query_point);        
         
         if( distance < min_distance and distance > FT(0) )
         {
             min_distance = distance;
             p1 = query_point;
             p0 = mesh.point(vit);
         }
      }
      return std::make_pair<Point_3&,Point_3&>(p0,p1);
   }

  // DocString: span
  /**
   * @brief Finds the span in a given Cartesian direction, i.e. (min, max) vertex position
   * 
   * @param direction Cartesian direction ( 0=x ,1=y,2=z) 
   * @returns a pair of doubles like (min,max)
   */
   std::pair<double,double> span(int direction)
   {
     assert_non_empty_mesh();
     auto bbox_3 = CGAL::Polygon_mesh_processing::bbox(mesh);
     std::pair<double,double> span (bbox_3.min(direction),bbox_3.max(direction));
     return span; 
   }

  /** 
   * @brief Adjust the vertex coordinates of vertices in a vector that it iterated over.
   * 
   * @tparam Inputiterator 
   * @param begin call std::vector
   * @param end call of std::vector
   * @param mulitplier of the normal vector that determines the vertes movement.
   */
   template<typename InputIterator >
   void adjust_vertices_in_region(InputIterator begin , InputIterator  end, const double c)
   { 
      assert_non_empty_mesh();
      std::vector<std::pair<vertex_descriptor, Point_3> > adjust; 
      for(; begin != end; ++begin)
      {
         Vector_3 delta= CGAL::Polygon_mesh_processing::compute_vertex_normal(*begin,mesh); 
         Point_3 p = mesh.point(*begin) + c*delta;
         adjust.push_back(std::make_pair(*begin, p));
      }
     for(std::pair<vertex_descriptor, Point_3> s : adjust)
        mesh.point(s.first) = s.second;
   }

   template<typename InputIterator >
   void adjust_vertices_in_region(InputIterator begin ,InputIterator  end, const Vector_3 disp)
   { 
      assert_non_empty_mesh();
      std::vector<std::pair<vertex_descriptor, Point_3> > adjust; 
      for(; begin != end; ++begin)
      {
         Point_3 p = mesh.point(*begin) + disp;
         adjust.push_back(std::make_pair(*begin, p));
      }
     for(std::pair<vertex_descriptor, Point_3> s : adjust)
        mesh.point(s.first) = s.second;
   }   
   
  /**
   * @brief Smooths the vertices in the iterator input. This is done by taking the sum of the vector edges for each vertex, and
   * multiplied with the constant double input.  
   * 
   * @note removes all isolated vertices to avoid error in algorithm   
   * @tparam Inputerator 
   * @param template iterator begin.
   * @param template iterator end.
   * @param const double c 
   */
   template<typename InputIterator>
   void smooth_laplacian_region(InputIterator  begin , InputIterator  end ,const double c)
   {
      assert_non_empty_mesh();
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
      std::vector<std::pair<vertex_descriptor, Point_3> > smoothed;
      for(; begin != end; ++begin)
      {
         Point_3 current = mesh.point(*begin);
         Vector_3 delta=CGAL::NULL_VECTOR;
         HV_const_circulator vbegin(mesh.halfedge(*begin),mesh), done(vbegin);
         do
         {
            delta += Vector_3(current,mesh.point(*vbegin));
            *vbegin++;
         }while(vbegin!=done);
         Point_3 p = current + c*delta/mesh.degree(*begin); 
         smoothed.push_back(std::make_pair(*begin, p));
      }
      for(std::pair<vertex_descriptor, Point_3> s : smoothed)
         mesh.point(s.first) = s.second;
   }

  /**
   * @brief Smooths the vertices in the iterator input. This is done by taking the sum of the vector edges for each vertex, and
   * multiplied with the constant double input.  
   * 
   * @param vertex_vector_map iterator begin.
   * @param vertex_vector_map iterator end.
   * @param const double c 
   * @overload
   */
   void smooth_laplacian_region(vertex_vector_map::iterator begin, vertex_vector_map::iterator end ,const double c)
   {
      assert_non_empty_mesh();
      CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
      std::vector<std::pair<vertex_descriptor, Point_3> > smoothed;
      for(; begin != end; ++begin)
      {
         Point_3 current = mesh.point(begin->first);
         Vector_3 delta=CGAL::NULL_VECTOR;
         HV_const_circulator vbegin(mesh.halfedge(begin->first),mesh), done(vbegin);
         do
         {
            delta += Vector_3(current,mesh.point(*vbegin));
            *vbegin++;
         }while(vbegin!=done);

         Point_3 p = current + c*delta/mesh.degree(begin->first);  // WHY C
         smoothed.push_back(std::make_pair(begin->first, p));
      } 
      for(std::pair<vertex_descriptor, Point_3> s : smoothed)
         mesh.point(s.first) = s.second;
   }

  // DocString: smooth_taubin
  /** 
   * @brief Taubin smothing of surface mesh. 
   *
   * Taubin smoothing of the surface vertices. This corresponds to a 
   * Laplacian smoothing with value \lambda, followed by a Laplacian smoothing with value \mu
   * Given the requirement:     
   *               $ \lambda < -\mu $    
   * @note The Laplacian smoothing parameters are set, but the 
   *       user may construct with their own parameters with smooth_laplacian function.
   * @param nb_iter number of iterations of smoothing   
   */
   void smooth_taubin(const size_t nb_iter) 
   {
       for(size_t i = 0; i < nb_iter; ++i) 
       {
           this->smooth_laplacian(0.8,1);  
           this->smooth_laplacian(-0.805,1);
       }
   }

  /** 
   * @brief Taubin smothing of specified surface mesh vertices 
   *
   * Taubin smaoothing of the input vertices. This corresponds to a 
   * Laplacian smoothing with value \lambda, followed by a Laplacian 
   * smoothing with value \mu: Given the requriment :     
   *               $ \lambda < -\mu $    
   *
   * @note The Laplacian smoothin parameters are set, but the 
   *       user can combine smooth_laplacian. 
   * @see [Articel]( https://www.doi.org/10.1109/ICCV.1995.466848) 
   * @param template iterator begin 
   * @param template iterator end 
   * @param nb_iter number of iterations of smoothing  
   */
   template<typename InputIterator > 
   void smooth_taubin_region(InputIterator begin , InputIterator end ,const size_t nb_iter)
   {
      assert_non_empty_mesh();
      for(size_t i = 0; i < nb_iter; ++i) 
      {
          this->smooth_laplacian_region(begin,end,0.8);
          this->smooth_laplacian_region(begin,end,-0.805);
      }
   }

  // DocString: num_self_intersections
  /** 
   * @brief Returns number of self intersection in the surface 
   * @returns number number of self intersection in the surface.
   */ 
   int num_self_intersections() 
   {
      assert_non_empty_mesh();
      std::vector< std::pair<face_descriptor, face_descriptor> > intersected_tris;
      CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(intersected_tris));
      return intersected_tris.size(); 
   }

  /** 
   * @brief Returns self-intersection triangles in the surface. 
   * @returns intersecting a vector of self-intersection triangles in the surface.
   */ 
   face_vector get_self_intersecting_faces() 
   {
      assert_non_empty_mesh();
      std::vector< std::pair<face_descriptor, face_descriptor> > intersected_tris;
      CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(intersected_tris));
      face_set collector; 
      for( auto tri : intersected_tris) 
      {
         collector.insert(tri.first);
         collector.insert(tri.second);
      }  
      face_vector intersecting(collector.begin(),collector.end());  
      return intersecting; 
   }

  /** 
   * @brief Returns vertices of triangles.
   * @param faces a vector of triangles . 
   * @returns vertex_vector associated faces. 
   */    
   vertex_vector get_vertices( face_vector faces) 
   {
      vertex_set collector;
      for( auto fd : faces ) 
      {
         CGAL::Vertex_around_face_circulator<Mesh> vbegin(mesh.halfedge(fd), mesh), done(vbegin);
         do
         {
            collector.insert(*vbegin);
                *vbegin++;
         }while(vbegin!=done);  
      }
      vertex_vector output(collector.begin(),collector.end()); 
      return output; 
   }

  /** 
   * @brief Returns edges of triangles.
   * @param faces a vector of triangles . 
   * @returns std::vector<edge_descriptor>  containing edges of triangles. 
   */    
   std::vector<edge_descriptor> get_edges(face_vector faces)
   {
      std::set< edge_descriptor> collector;
      for( auto fd : faces ) 
      {
         CGAL::Halfedge_around_target_circulator<Mesh> hbegin(mesh.halfedge(fd),mesh), done(hbegin);
         do
         {
            collector.insert(mesh.edge(*hbegin));
            *hbegin++;
         }while(hbegin!=done);  
      
      }
      std::vector<edge_descriptor> output(collector.begin(),collector.end());   
      return output;
   }
   
  /** 
   * @brief Removes faces from surface mesh with CGAL::Euler.
   * @param faces a vector of triangles . 
   * @returns none
   */      
   void remove_faces(face_vector faces) 
   {
       for( auto fd : faces) 
            CGAL::Euler::remove_face(mesh.halfedge(fd),mesh);
        this->fill_holes();
   }
   
  // DocString: collapse_edges
  /**
   * @brief Combines smaller edges together, so that all edges are larger than the input parameter
   * @see [edge_collapse](https://doc.cgal.org/latest/Surface_mesh_simplification/index.html)
   * @param  target_edge_length the lower bound of the edges in the surface mesh.
   * @returns r the number of collapsed edges 
   */
   int collapse_edges(const double target_edge_length)
   {
      assert_non_empty_mesh();
      CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<double> stop(target_edge_length);

      const int r = CGAL::Surface_mesh_simplification::edge_collapse(
                    mesh,
                    stop,
                    CGAL::parameters::get_cost(CGAL::Surface_mesh_simplification::Edge_length_cost<Mesh>())
                    .get_placement(CGAL::Surface_mesh_simplification::Midpoint_placement<Mesh>()));
      return r;
   }

  // DocString: collapse_edges
  /**
   * @brief Collapses smaller edges together, given that the edges have the same 
   *        direction.
   * @returns the number of collapsed edges
   * @overload  
   */
   int collapse_edges() 
   {
      Cost_stop_predicate<Mesh> stop(1.e-6);
      const int r = CGAL::Surface_mesh_simplification::edge_collapse(mesh,stop);
      return r;
   }

  // DocString: get_slice  
  /**
   * @brief Slices a SVMTK Surface object according to a given plane, 
   *        and returns a SVMTK Slice object.
   *
   * @see [Plane_3] (https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html)  
   * @see [Slice](Slice.h) 
   * @tparam SVMTK Slice class. 
   * @param plane_3 a wrapped CGAL Plane_3 object 
   * @returns slice SVMTK Slice object.
   */
   template<typename Slice>
   std::shared_ptr<Slice> get_slice(Plane_3 plane_3)  
   {
      assert_non_empty_mesh();  

      CGAL::Polygon_mesh_slicer<Mesh, Kernel> slicer(mesh); 
      Polylines polylines_3D;
      slicer(plane_3, std::back_inserter(polylines_3D));

      std::vector<std::vector<Point_2>> polylines_2;
      for(auto pol = polylines_3D.begin(); pol != polylines_3D.end(); ++pol) 
      {
         std::vector<Point_2> result;
         std::vector<Point_2> polyline_2;
         for(auto pit = pol->begin(); pit != pol->end(); ++pit)
             polyline_2.push_back(plane_3.to_2d(*pit));
         polylines_2.push_back(polyline_2);
      }
      std::shared_ptr<Slice> slice(new Slice(plane_3,polylines_2)); 
      return slice;
   }  

  // DocString: get_slice
  /**
   * @brief  Slices a surface mesh based on a plane, that is defined by the plane equation : 
   *       p1*x1 + p2*x2 +p3*x3 -x4 = 0.
   *
   * @note Better result if the plane normal is close to unity. 
   * @param x1 plane equation parameter.
   * @param x2 plane equation parameter.
   * @param x3 plane equation parameter.
   * @param x4 plane equation parameter.
   * @returns a SVMTK slice object. 
   * @overload
   */
   template<typename Slice>
   std::shared_ptr<Slice> get_slice(double x1,double x2, double x3 ,double x4)
   {
      assert_non_empty_mesh();
      if( x1==0 && x2==0 && x3==0 )
         throw InvalidArgumentError("Invalid plane parameters.");
      Plane_3 plane = Plane_3(x1, x2, x3, x4);
      return this->get_slice<Slice>(plane);
   }

  /**
   * @brief  Slices a SVMTK Surface object according to a plane, 
   * and returns the intersecting lines as a vetor of point vectors
   * @see [Plane_3] (https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html)   
   * @returns a SVMTK Slice object.
   */
   Polylines polylines_in_plane(Plane_3 plane_3) 
   {
      assert_non_empty_mesh();
      CGAL::Polygon_mesh_slicer<Mesh, Kernel> slicer(mesh); 
      Polylines polylines;
      slicer(plane_3, std::back_inserter(polylines));
      return polylines;
   } 
     
  /**
   * @brief  Slices a SVMTK Surface object according to a plane.
   * and returns the intersecting lines as a vetor of point vectors.
   * @see [Plane_3] (https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html)   
   * @returns a SVMTK Slice object.
   * @overload
   */   
   Polylines polylines_in_plane(Point_3 p, Vector_3 v)
   {
      return polylines_in_plane(Plane_3(p,v));
   } 
  
  // DocString: isotropic_remeshing   
  /** 
   * @brief Isotropic remeshing of surface mesh. Remeshing of the surface mesh so that all edges have the same length.
   * 
   * Uses CGAL isotropic_remeshing and split_long_edges. 
   * @note The option protect_border, then edges over the threshold of 40 degrees is found and protected. 
   * @see [split_long_edges](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html) 
   * @see [isotropic_remeshing](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html) 
   * @note split_long_edges is used to avoid a pitfall described [here](https://doc.cgal.org/5.0.3/Polygon_mesh_processing/index.html)
   *
   * @param target_edge_length the edge length that is targeted in the remeshed patch. If 0 is passed then only the edge-flip, tangential relaxation, and projection steps will be done. 
   * @param nb_iter the number of iterations for the sequence of atomic operations performed. 
   * @param protect_border If true, constraint edges cannot be modified at all during the remeshing process. 
   */
   void isotropic_remeshing(double target_edge_length, unsigned int nb_iter, bool protect_border)
   {
       assert_non_empty_mesh();
       
       if( protect_border )
       {  
          EIFMap eif = get(CGAL::edge_is_feature, mesh);
          PIMap pid = get(CGAL::face_patch_id_t<int>(), mesh);

          CGAL::Polygon_mesh_processing::sharp_edges_segmentation(mesh, 60, eif, pid);
       
          CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length, mesh, CGAL::parameters::edge_is_constrained_map(eif));
          CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh),
                              target_edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::edge_is_constrained_map(eif)
                              .number_of_iterations(nb_iter)
                              .protect_constraints(protect_border));      

       }
       else 
       {
       CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
       CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh),
                              target_edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter));
       }                         
                              
   }

   void isotropic_remeshing(face_vector faces, double target_edge_length, unsigned int nb_iter, bool protect_border)
   {
       assert_non_empty_mesh();
       CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
       CGAL::Polygon_mesh_processing::isotropic_remeshing(faces,
                              target_edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
                              .protect_constraints(protect_border));
   }

  // DocString: clip
  /**
   * @brief Clips the surface mesh.
   * Clips a surface mesh based on a plane defined with the equation 
   *                    $ ax+by+cz +d = 0  $
   * @see [CGAL::Polygon_mesh_processing::clip](https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#Coref_section)
   * @see [CGAL::Plane_3](https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html)
   * @note Setting plane equation parameters close to unit values produce better clip.
   *
   * @param  a parameter in plane equation.
   * @param  b parameter in plane equation.
   * @param  c parameter in plane equation.
   * @param  d parameter in plane equation.
   * @param  preserve_manifold true to preserve manifold.  
   * @returns true if successful.
   */
   bool clip(double a,double b, double c ,double d, bool preserve_manifold)
   {      
     assert_non_empty_mesh();
     return clip(Plane_3(a,b,c,d), preserve_manifold);
   }

   
  // DocString: clip
  /**
   * @brief Clips a surface mesh given CGAL Plane_3 object.
   *
   * @param  plane Wrapped CGAL Plane_3 object.
   * @param  preserve_manifold true to preserve manifold.  
   * @overload 
   */
   bool clip(Plane_3 plane, bool preserve_manifold)
   {
      assert_non_empty_mesh();
      Point_3 center = this->centeroid();
      Point_3 point = plane.projection(center);
      Vector_3 vector = plane.orthogonal_vector();
      return clip(point, vector, preserve_manifold);
   }

  // DocString: clip
  /**
   * @brief Clips a surface mesh given another SVMTK Surface object. 
   * @param  other SVMTK Surface object.
   * @param  preserve_manifold true to preserve manifold. 
   * @param invert if true invert other surface.  
   * @overload
   */
  bool clip(Surface other, bool invert, bool preserve_manifold)
  {
    assert_non_empty_mesh();
    if(invert) 
       CGAL::Polygon_mesh_processing::reverse_face_orientations(other.mesh);
       
    double target_edge_length = average_edge_length();  
    if( other.does_bound_a_volume() ) 
       return CGAL::Polygon_mesh_processing::clip(mesh, other.mesh,CGAL::Polygon_mesh_processing::parameters::clip_volume(preserve_manifold) );

    auto sucess = CGAL::Polygon_mesh_processing::clip(mesh, other.mesh, CGAL::Polygon_mesh_processing::parameters::clip_volume(false));

    if( preserve_manifold )
    {
        isotropic_remeshing(target_edge_length,5,false);

        std::vector<halfedge_descriptor> border_cycles;
        CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
      
        std::vector<face_descriptor> patch;
        for( halfedge_descriptor h : border_cycles ) 
             CGAL::Polygon_mesh_processing::triangulate_hole(mesh, h, CGAL::parameters::face_output_iterator(std::back_inserter(patch)));
        
        CGAL::Polygon_mesh_processing::isotropic_remeshing(patch,
                              target_edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));     
    }
    return sucess;	
  } 
  
  // DocString: clip 
  /**
   * @brief Clips the surface mesh with circle in plane
   *
   * Creates a circle in a specified plane, and uses this 
   * circle to clip the surface mesh. 
   *  
   * @param  point on the plane 
   * @param  vector a vector normal to the plane.
   * @param  radius of the circle 
   * @param  invert the clip
   * @param  preserve_manifold true to preserve manifold  
   * @returns true if succesfull 
   */
   bool clip(Point_3 point, Vector_3 vector, double radius, bool invert, bool preserve_manifold)
   {
     Surface circle;
     circle.make_circle_in_plane(point, vector, radius, average_edge_length() ); 
     return clip(circle, invert, preserve_manifold);
   }
    
  // DocString: clip 
  /**
   * @brief Clips the surface mesh with a circle that has the minimum bounding radius, gives better result 
   *
   * Creates a circle in a specified plane, and uses this 
   * circle to clip the surface mesh. 
   *  
   * @param  point on the plane 
   * @param  vector a vector normal to the plane.
   * @param  preserve_manifold true to preserve manifold  
   * @returns true if succesfull 
   */                        
   bool clip(Point_3 point, Vector_3 vector, bool preserve_manifold)
   {   
     assert_non_empty_mesh();
     Surface clipper;
     double radius = get_bounding_radius(); 
     double edge_length = average_edge_length();
     clipper.make_circle_in_plane(point, vector,radius, edge_length); 
     return clip(clipper, true, preserve_manifold);
   }
   
  //DocString: get_perpendicular_cut
  /**
   * @brief Constructs a clip surface where a given plane intersects the mean curvature flow[], the clip surface will be  
   *        perpendicular to mean curvature flow.
   * @param plane that intersects the mean curvature flow (centerlines)  
   * @param radius of the clip surface.
   * @return SVMTK Surface object, circle in plane 
   */
   std::shared_ptr<Surface> get_perpendicular_cut(Plane_3 plane, double radius)
   {
        Polylines polylines = this->mean_curvature_flow();
        
        bool _do_intersect = false;
        Vector_3 normal; 
        Point_3 intersecting_point;
        
        for ( auto line : polylines   ) 
        {
            // Make function
            for ( auto pit1 = line.begin(), pit2 = std::next(line.begin()); pit2 != line.end(); pit1++, pit2++)
            {
               if ( plane.has_on_positive_side(*pit1) != plane.has_on_positive_side(*pit2) )
               {
                    intersecting_point = plane.projection(*pit2);
                    normal = Vector_3(*pit1,*pit2);
                    normal/=CGAL::sqrt(normal.squared_length());
                     _do_intersect = true;
                    const auto result = CGAL::intersection(Line_3(*pit1,*pit2), plane);
                    if (result) 
                    {
                       const Point_3* p = std::get_if<Point_3 >(&*result);
                       intersecting_point = *p; 
                    } 
                    break;
               }             
            }
            if ( _do_intersect ) 
               break;
        }
        
        if ( !_do_intersect )   
             throw InvalidArgumentError("Surface and plane does not intersect for a perpendicular cut.");    
        else 
        {
        
        std::shared_ptr<Surface> result(new Surface()); 
        double ael = average_edge_length();
        
        if ( radius == 0.0 ) 
            radius = distance_to_point( intersecting_point) + 4.0 ;
        
        // Principal direction should be positive for consistiency
        if ( normal.x() + normal.y() + normal.z() < 0)
            normal = -normal;
        
        result->make_circle_in_plane(intersecting_point, normal,radius, ael);    
        return result;
        }
       
   }
   
  //DocString: get_perpendicular_cut
  /**
   * @brief Constructs a clip surface near a query point. The clip surface will be  
   *        perpendicular to mean curvature flow.
   *
   * @param query SVMTK Point_3 object   
   * @param radius of the clip surface.
   * @return SVMTK Surface object, circle in plane 
   */ 

   std::shared_ptr<Surface> get_perpendicular_cut(Point_3 query, double  radius)
   {
        Polylines polylines = this->mean_curvature_flow();
        Vector_3 normal; 
        Point_3 closest;
        FT min_dist(100);
        FT sq_dist; 
        
        for ( auto line : polylines   ) 
        {
                for ( auto pit1 = line.begin(), pit2 = std::next(line.begin()); pit2 != line.end(); pit1++, pit2++)
                {
                    sq_dist = CGAL::squared_distance(query,*pit1);
                                                         
                    if (sq_dist <  min_dist ) 
                    {
                        min_dist = sq_dist;
                        closest = *pit1;
                        normal = Vector_3(*pit1,*pit2);
                    }
                }

                sq_dist = CGAL::squared_distance(query,*line.end() );
                if (sq_dist <  min_dist ) 
                {
                        min_dist = sq_dist;
                        closest = *line.end() ;
                        normal = Vector_3( *line.end() ,*std::prev(line.end()));
                }
                
        } 
        if ( min_dist > FT(100) ) 
           throw InvalidArgumentError("Point is not close enough to surface.");      
        std::shared_ptr<Surface> result(new Surface()); 
        double ael = average_edge_length();
        if ( radius == 0.0 ) 
             radius = distance_to_point( closest) + 4.0 ;
        result->make_circle_in_plane(closest, normal, radius, ael);    
        return result;        
   }

  // DocString: fill_holes
  /**
   * @brief  Finds and fills holes in surface mesh.
   *
   * Uses CGAL function triangulate_refine_and_fair_hole. 
   * @see [triangulate_refine_and_fair_hole](https://doc.cgal.org/latest/Polygon_mesh_processing/group__hole__filling__grp.html)
   * @returns the nb_holes number of holes filled.
   */
   std::pair<bool,int> fill_holes()
   {
     unsigned int nb_holes = 0;
     CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
     CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices(mesh);
     
     std::vector<halfedge_descriptor> border_cycles;
     
     CGAL::Polygon_mesh_processing::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));

     for(halfedge_descriptor h : border_cycles)
     {
         std::vector<face_descriptor>  patch_facets;
         std::vector<vertex_descriptor> patch_vertices;
         bool success = std::get<0>(CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(mesh, h,
                                                 CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets)).
                                                                   vertex_output_iterator(std::back_inserter(patch_vertices))));
         if( success )                                                               
           nb_holes++;

     }
     set_outward_face_orientation(); 
     return std::make_pair(does_bound_a_volume(),nb_holes);
   }

  // DocString: triangulate_faces
  /** 
   * @brief  Triangulates faces of the surface mesh.
   * Uses CGAL function triangulate_faces.
   * @see [triangulate_faces] (https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html)
   * @returns true if complete.
   */
   bool triangulate_faces()
   {
      assert_non_empty_mesh();
      CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
      BOOST_FOREACH(face_descriptor fd, faces(mesh))
           if( next(next(halfedge(fd, mesh), mesh), mesh)!=prev(halfedge(fd, mesh), mesh) )
              std::cerr << "Error: non-triangular face left in mesh." << std::endl;
     return true;
   }

  /**
   * @brief Returns the points corresponding to vertices of the surfaces mesh.
   * @param vertices a vector of vertices. 
   * @returns a vector of surface points.
   */
   point_vector get_points(vertex_vector &vertices) 
   {
      assert_non_empty_mesh();
      point_vector result;
      for(vertex_descriptor vit : vertices)
          result.push_back(mesh.point(vit));
      return result;
   } 

   Point_3 get_point(vertex_descriptor vd) 
   {
      assert_non_empty_mesh();
      return mesh.point(vd);
   } 
   
  // DocString: get_points
  /**
   * @brief Returns the points of the surfaces mesh.
   * @returns vector of surface points.
   */ 
   point_vector get_points() 
   {
     assert_non_empty_mesh();
     point_vector result;
     for(vertex_descriptor vit : mesh.vertices())
          result.push_back(mesh.point(vit));
     return result;
   }  

  // DocString: is_point_inside
  /**
   * @brief Checks if a point is inside surface mesh.
   *  
   * Query if a point is inside the member surface mesh.
   * @param point_3 surface mesh point.
   * @returns bool true if point is inside surface otherwise false. 
   */
   bool is_point_inside(Point_3 point_3) 
   {
      assert_non_empty_mesh();
      Inside is_inside_query(mesh);
      CGAL::Bounded_side res = is_inside_query(point_3);
      if (res == CGAL::ON_BOUNDED_SIDE or res == CGAL::ON_BOUNDARY)
         return true;
      else 
         return false;
   }
   
   bool is_point_on_boundary(Point_3 point_3) 
   {
      assert_non_empty_mesh();
      Inside is_inside_query(mesh);
      CGAL::Bounded_side res = is_inside_query(point_3);
      if ( res==CGAL::ON_BOUNDARY )
         return true;
      else 
         return false;
   }   

  // DocString: adjust_boundary
  /**
   * @brief Adjust surface mesh vertex coordinates in the normal vertex direction multiplied with an 
   * argument value.
   * 
   * @param c double value multiplier.  
   */ 
   void adjust_boundary(const double c)
   {
      assert_non_empty_mesh();
      Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
      adjust_vertices_in_region(vb, ve,c);
   }
   
  // DocString: adjust_boundary
  /**
   * @brief Adjust surface mesh vertex coordinates in the normal vertex direction multiplied with an 
   * argument value.
   * 
   * @param c double value multiplier.  
   */    
   void adjust_boundary(Surface other,const double c)
   {
      assert_non_empty_mesh();
      vertex_vector vertices = get_vertices_with_property<CGAL::ON_BOUNDED_SIDE,CGAL::ON_BOUNDARY>(other); 

      adjust_vertices_in_region(vertices.begin(), vertices.end(), c);
   }

  // DocString: smooth_laplacian
  /** 
   * @brief Smooths all vertices of surface mesh with Laplacian smoothing.
   *  
   * This is done by taking the sum of the vector edges for each vertex, and
   * mmultiplied with the constant double input.  
   *
   * @param c a double multipler of the displacment vector that decides the new vertex coordinate.   
   * @param iter integer that decides the number of iterations of algorithm
   * @note zero and negative integer will have no effect, and should therefore be avoided.     
   *  
   */
   void smooth_laplacian(const double c, int nb_iter)
   {
      assert_non_empty_mesh();
      Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
      for(int i=0; i<nb_iter; ++i)
          this->smooth_laplacian_region(vb, ve,c);
   }

  // DocString: smooth_shape
  /**
   * @brief Smooth the triangulated surface. 
   * @see [CGAL::Polygon_mesh_processing::smooth_shape](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html)
   *
   * @param a time step that corresponds to the speed by which the surface is smoothed. A larger time step results in faster 
   *       convergence but details may be distorted to have a larger extent compared to more iterations with a smaller step. 
   *       Typical values scale in the interval (1e-6, 1].
   * @param nb_iter number of iteations.  
   */
   void smooth_shape(double time,int nb_iter)
   {
      assert_non_empty_mesh();
      CGAL::Polygon_mesh_processing::smooth_shape(mesh, time, CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter));
   }


     //move to utility
   struct sort_vectors_by_length {
                  template<typename T>
                  bool operator()(const std::vector<T> & a, const std::vector<T> & b)
                       { return ( length_polyline_3(a.begin(), a.end())>length_polyline_3(b.begin(), b.end()) );}                   
       };   
  // DocString: mean_curvature_flow 
  /** 
   * @brief Computes the centerline of the surface.
   * @see [CGAL::extract_mean_curvature_flow_skeleton](https://doc.cgal.org/latest/Surface_mesh_skeletonization/index.html)
   * @returns Polylines vector of points. 
   */
   Polylines mean_curvature_flow() 
   {
      assert_non_empty_mesh();
    
      struct Construct_polylines
      {   
         const Skeleton& skeleton;
         std::vector<Point_3> polyline_3;
         std::vector<std::vector<Point_3>> polylines_3;
 
         Construct_polylines(const Skeleton& skeleton): skeleton(skeleton){}
         void start_new_polyline()
         {
            polyline_3.clear();
         } 
         void add_node(Skeleton_vertex v)
         {
             polyline_3.push_back(skeleton[v].point);
         }
         void end_polyline()
         {
            polylines_3.push_back(polyline_3);                     
         }
         std::vector<std::vector<Point_3>>& get_polylines()
         {
             return polylines_3;
         }         
      };
  
      Polylines polylines_3;
      Skeleton skeleton;
      CGAL::extract_mean_curvature_flow_skeleton(mesh, skeleton);
      Construct_polylines Visitor(skeleton);
      CGAL::split_graph_into_polylines(skeleton,Visitor);
      auto lines = Visitor.get_polylines();
      
 
      
      std::sort(lines.begin(), lines.end(), sort_vectors_by_length());
      return lines;
  }
  
  // 
  /**
   *
   */
  Polylines mean_curvature_flow( std::string ordering) 
  {
      Polylines output = mean_curvature_flow();
      output = ordering_polylines_3(output,ordering);
      return output;
  
  }
  
  
  

  // DocString: get_shortest_surface_path
  /**
   * @brief Computes and returns the shortest surface path between two points.
   * Two points are projected on to the surface mesh, and the shortest surface path
   * is found.
   * @param source first point 
   * @param target second point.
   * @returns a vector of sequential points, i.e. polyline. 
   */
   Polyline get_shortest_surface_path(Point_3 source, Point_3 target) 
   {     
      assert_non_empty_mesh();

      std::vector<Point_3> points;

      AABB_Tree tree(faces(mesh).first, faces(mesh).second, mesh); 
      Surface_mesh_shortest_path shortest_paths(mesh);     

      Face_location source_location = shortest_paths.locate(source,tree);
      Face_location target_location = shortest_paths.locate(target,tree);

      shortest_paths.add_source_point(source_location);
      shortest_paths.shortest_path_points_to_source_points(target_location.first,target_location.second, std::back_inserter(points));  
     
      return points;
   } 
   
  // DocString: get_shortest_surface_path
  /**
   * @brief Computes and returns the shortest surface path between two points. 
   * @param x0 first point x coordinate
   * @param y0 first point y coordinate
   * @param z0 first point z coordinate
   * @param x1 second point x coordinate
   * @param y1 second point y coordinate
   * @param z1 second point z coordinate
   * @returns a polyline of the shortest surface path.
   * @overload
   */
   Polyline get_shortest_surface_path(double x0, double y0, double z0, double x1, double y1, double z1)
   {  
     Point_3 source(x0,y0,z0);
     Point_3 target(x1,y1,z1);
     return get_shortest_surface_path(source, target); 
   }

  // DocString: make_cube   
  /**
   * @brief Creates a surface mesh structure with vertices and facets connecting vertices for a cube.
   * 
   * @param x0 first point x coordinate
   * @param y0 first point y coordinate
   * @param z0 first point z coordinate
   * @param x1 second point x coordinate
   * @param y1 second point y coordinate
   * @param z1 second point z coordinate
   * @param edge_length the size of an edge in each direction
   */
   void make_cube( double x0, double y0, double  z0,  double x1, double y1, double z1, double edge_length) 
   {
      typedef boost::multi_array<int, 3> array_type;
     
      std::vector<face_vector> sides;
      clear(); 
      if( x0==x1 or y0==y1 or z0==z1 ) 
         throw InvalidArgumentError("Invalid argument.");
      double dx =(x1-x0);  
      double dy =(y1-y0);
      double dz =(z1-z0);
     
      //if( abs(dx)<edge_length || abs(dy)<edge_length ||  abs(dz)<edge_length) 
      //   throw InvalidArgumentError("Select smaller edge length.");    
      if( abs(dx)<edge_length || abs(dy)<edge_length ||  abs(dz)<edge_length) 
          edge_length = std::min( abs(dx), std::min(abs(dy), abs(dz)));


      int Nx = abs(static_cast<int>(dx/edge_length));
      int Ny = abs(static_cast<int>(dy/edge_length));
      int Nz = abs(static_cast<int>(dz/edge_length));
      double ddx = dx/static_cast<double>(Nx); 
      double ddy = dy/static_cast<double>(Ny); 
      double ddz = dz/static_cast<double>(Nz);  
      
      array_type map(boost::extents[Nx+1][Ny+1][Nz+1]);
         
      int index =0;
      for(int i=0; i<Nx+1; ++i)
      {
          for(int j=0; j<Ny+1; ++j)
          { 
              for(int k=0; k<Nz+1; ++k)
              {
                  if( i==0 or i==Nx or j==0 or j==Ny or k==0 or k==Nz )
                  {     
                     mesh.add_vertex(
                          Point_3(x0+static_cast<double>(i)*ddx,
                                  y0+static_cast<double>(j)*ddy,
                                  z0+static_cast<double>(k)*ddz) 
                     );
                     map[i][j][k] = index++;    
                  }     
              }
         }    
      }
      face_vector s1;
      face_vector s2;
      for(int j=0; j<Ny; ++j)
      { 
          for(int k=0; k<Nz; ++k)
          {
              s1.push_back(mesh.add_face(Index(map[0][j][k]),    Index(map[0][j+1][k+1]),  Index(map[0][j+1][k])));
              s1.push_back(mesh.add_face(Index(map[0][j][k]),    Index(map[0][j][k+1]),    Index(map[0][j+1][k+1])));
              s2.push_back(mesh.add_face(Index(map[Nx][j][k]),   Index(map[Nx][j+1][k]),   Index(map[Nx][j][k+1])));
              s2.push_back(mesh.add_face(Index(map[Nx][j+1][k]), Index(map[Nx][j+1][k+1]), Index(map[Nx][j][k+1])));
          }
      }
      face_vector s3;
      face_vector s4;
      for(int i=0; i<Nx; ++i)
      { 
          for(int k=0; k<Nz; ++k)
          {
              s3.push_back(mesh.add_face(Index(map[i][0][k]),    Index(map[i+1][0][k]),    Index( map[i][0][k+1])));
              s3.push_back(mesh.add_face(Index(map[i+1][0][k]),  Index( map[i+1][0][k+1]), Index(map[i][0][k+1])));
              s4.push_back(mesh.add_face(Index(map[i][Ny][k]),   Index(map[i][Ny][k+1]),   Index(map[i+1][Ny][k+1])));
              s4.push_back(mesh.add_face(Index(map[i+1][Ny][k]), Index(map[i][Ny][k]),     Index(map[i+1][Ny][k+1])));
          } 
      }
      face_vector s5;
      face_vector s6;
      for(int i=0; i<Nx; ++i)
      { 
          for(int j=0; j<Ny; ++j)
          {
              s5.push_back(mesh.add_face(Index(map[i][j][0]),  Index(map[i][j+1][0]),    Index(map[i+1][j+1][0])));
              s5.push_back(mesh.add_face(Index(map[i][j][0]),  Index(map[i+1][j+1][0]),  Index(map[i+1][j][0])));
              s6.push_back(mesh.add_face(Index(map[i][j][Nz]), Index(map[i+1][j][Nz]),   Index(map[i+1][j+1][Nz])));
              s6.push_back(mesh.add_face(Index(map[i][j][Nz]), Index(map[i+1][j+1][Nz]), Index(map[i][j+1][Nz])));
          } 
      }    
      sides.push_back(s1);
      sides.push_back(s2);
      sides.push_back(s3);
      sides.push_back(s4);
      sides.push_back(s5);
      sides.push_back(s6);     
     
      set_outward_face_orientation(); 

      for(auto side : sides) 
      {     
          CGAL::Polygon_mesh_processing::isotropic_remeshing(side,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));        
      }
      return;
   }
   
  // DocString: make_cube   
   /**
    *  @brief Constructs a cube surface.
    *  @param p0 SVMTK 3D point of a corner in a cube 
    *  @param p1 SVMTK 3D point of the corner opposite of p0.
    *  @param edge_length the upper edge length of the constructed triangulated surface 
    *  @overload 
    */ 
   void make_cube(Point_3 p0, Point_3 p1, double edge_length )
   {
        make_cube(CGAL::to_double(p0.x()), CGAL::to_double(p0.y()),CGAL::to_double(p0.z()),
                   CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),CGAL::to_double(p1.z()),
                   edge_length);
   }
   
  // DocString: make_cylinder
  /**
   * @brief Creates a surface mesh structure with vertices and facets connecting vertices for a cylinder
   * @param x0 x coordinate of the first cylinder center 
   * @param y0 y coordinate of the first cylinder center
   * @param z0 z coordinate of the first cylinder center
   * @param x1 x coordinate of the second cylinder center 
   * @param y1 y coordinate of the second cylinder center
   * @param z1 z coordinate of the second cylinder center
   * @param r0 the radius of the cylinder 
   * @param edge_length the target length of an edge in the mesh.
   */
   void make_cylinder(double x0, double y0, double z0, double x1, double y1, double z1, double r0, double edge_length)
   {    
      clear();
      this->make_cone(x0, y0, z0, x1, y1, z1, r0, r0, edge_length);
   }
   
   // DocString: make_cylinder   
   /**
   *  @brief Creates a surface mesh structure with vertices and facets connecting vertices for a cylinder.
   *  @param p0 SVMTK 3D point of the bottom center of the cylinder.
   *  @param p1 SVMTK 3D point of the top center of the cylinder.
   *  @param r0 the radius of the cylinder. 
   *  @param edge_length the upper edge length of the constructed triangulated surface 
   *  @overload 
   */
   void make_cylinder(Point_3 p0,Point_3 p1, double r0, double edge_length)
   { 
        make_cylinder(CGAL::to_double(p0.x()), CGAL::to_double(p0.y()),CGAL::to_double(p0.z()),
                       CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),CGAL::to_double(p1.z()),
                       r0,edge_length);
   }

  /**
   * @brief Creates a surface mesh structure with vertices and facets connecting vertices for a cone with
   * radius equal zero.
   * 
   * @param x0 x coordinate of the first cylinder center. 
   * @param y0 y coordinate of the first cylinder center.
   * @param z0 z coordinate of the first  cylinder center.
   * @param x1 x coordinate of the second cylinder center.
   * @param y1 y coordinate of the second cylinder center.
   * @param z1 z coordinate of the second cylinder center.
   * @param r0 the none-zero radius of the cone, corresponding to the first cylinder center.  
   * @param edge_length the target length of an edge in the mesh.
   */
   void make_cone_( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0, double edge_length) 
   {
      clear();
      face_vector fv1,fv3;
      Index v0 = mesh.add_vertex(Point_3(x0,y0,z0));
      Index v1 = mesh.add_vertex(Point_3(x1,y1,z1));        
      Index vb, vt;  
  
      Vector_3 dir(Point_3(x0,y0,z0), Point_3(x1,y1,z1));
      double l1 = std::sqrt(dir.squared_length());
      Plane_3 plane(Point_3(x1,y1,z1), dir/l1); 
    
      Vector_3 t1 = plane.base1();
      Vector_3 t2 = plane.base2();    
      t1=t1/std::sqrt(t1.squared_length()); 
      t2=t2/std::sqrt(t2.squared_length()); 
     
      int number_of_segments = static_cast<int>(CGAL_PI*2*r0/edge_length);    
    
      if( number_of_segments<3) 
         throw InvalidArgumentError("Select smaller edge length.");    
        
      double c = 360.0/(double)number_of_segments;
      for(int i = 0; i<number_of_segments; ++i) 
      {              
          Point_3 pb= Point_3(x0,y0,z0)  + t1*r0*std::cos(c*i*CGAL_PI/180) + t2*r0*std::sin(c*i*CGAL_PI/180); 
          vb = mesh.add_vertex(pb);
          if( i!=0 )
          {  
            fv1.push_back( mesh.add_face(v0, vb, Index(vb-1)));
            fv3.push_back( mesh.add_face(Index(vb-1), vb, v1));
          }
      }
      fv1.push_back(mesh.add_face(Index(2), vb, v0));
      fv3.push_back(mesh.add_face(vb, Index(2), v1));

      CGAL::Polygon_mesh_processing::isotropic_remeshing(fv1,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));

      CGAL::Polygon_mesh_processing::isotropic_remeshing(fv3,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));
   }

  // DocString: make_cone
  /**
   * @brief Creates a surface mesh structure with vertices and facets connecting vertices for a cone.
   *
   * @note The function also handles the special cases of sharp cone and cylinder. 
   * @param x0 x coordinate of the first cone center 
   * @param y0 y coordinate of the first cone  center
   * @param z0 z coordinate of the first cone center
   * @param x1 x coordinate of the second cone center 
   * @param y1 y coordinate of the second cone center
   * @param z1 z coordinate of the second cone center
   * @param r0 the radius corresponding to the first cone center.   
   * @param r1 the radius corresponding to the second cone center. 
   * @param edge_length the target length of an edge in the mesh.
   */
   void make_cone(double x0, double y0, double  z0, double x1, double y1, double z1, double r0, double r1, double edge_length) 
   {  
      if( (r0*r0 +r1*r1)==0 ) 
         throw InvalidArgumentError("Select one non-zero radius.");   
         
      clear();
      if( r0==0.0 )
      {
         make_cone_(x1,y1,z1,x0,y0,z0,r1,edge_length);
         return;
      }
      else if( r1==0.0 )
      {
         make_cone_(x0,y0,z0,x1,y1,z1,r0,edge_length);
         return;
      }
      face_vector fv1,fv2,fv3;

      Index v0 = mesh.add_vertex(Point_3(x0,y0,z0));
      Index v1 = mesh.add_vertex(Point_3(x1,y1,z1));      
      Index vb, vt;
      Vector_3 dir(Point_3(x0,y0,z0),Point_3(x1,y1,z1));
      double l1 = std::sqrt(dir.squared_length());
      
      Plane_3 plane(Point_3(x1,y1,z1),dir/l1);
      
      Vector_3 t1 = plane.base1();
      Vector_3 t2 = plane.base2();
      
      t1=t1/std::sqrt(t1.squared_length());    
      t2=t2/std::sqrt(t2.squared_length()); 
     
      double r = std::max(r0,r1);
     
      int number_of_segments = static_cast<int>(CGAL_PI*2*r/edge_length); // approx pi, since setting it as int    
     
      if( number_of_segments<3 ) 
         throw InvalidArgumentError("Select smaller edge length.");    
     
      double c  = 360.0/(double)number_of_segments;
     
      for(int i = 0; i < number_of_segments; ++i) 
      {         
          Point_3 pb = Point_3(x0,y0,z0)  + t1*r0*std::cos(c*i*CGAL_PI/180) + t2*r0*std::sin(c*i*CGAL_PI/180); 
          Point_3 pt = Point_3(x1,y1,z1)  + t1*r1*std::cos(c*i*CGAL_PI/180) + t2*r1*std::sin(c*i*CGAL_PI/180); 
          vb = mesh.add_vertex(pb);
          vt = mesh.add_vertex(pt);
          if( i!=0 )
          {  
            fv1.push_back( mesh.add_face(v0, vb,Index(vb-2)));
            fv2.push_back( mesh.add_face(v1, Index(vt-2), vt));
            fv3.push_back( mesh.add_face(Index(vt-2), Index(vb-2), vt));
            fv3.push_back( mesh.add_face(vb, vt, Index(vb-2)));
          }
      }
      fv1.push_back(mesh.add_face(Index(0), Index(2), vb));
      fv2.push_back(mesh.add_face(Index(1), vt, Index(3)));
      fv3.push_back(mesh.add_face(vt, vb, Index(3)));
      fv3.push_back(mesh.add_face(Index(2), Index(3), vb));

      CGAL::Polygon_mesh_processing::isotropic_remeshing(fv1,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));

      CGAL::Polygon_mesh_processing::isotropic_remeshing(fv2,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)           
                              .protect_constraints(true));
                              
      CGAL::Polygon_mesh_processing::isotropic_remeshing(fv3,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));
   }

  // DocString: make_cone  
  /**
   *  @brief Constructs a cone surface.
   *  @param p0 SVMTK 3D point of the bottom center of the cone
   *  @param p1 SVMTK 3D point of the top center of the cone.
   *  @param r0 the bottom radius of the cone. 
   *  @param r1 the top radius of the cone.     
   *  @param edge_length the upper edge length of the constructed triangulated surface 
   *  @overload 
   */
   void make_cone(Point_3 p0,Point_3 p1, double r0, double r1, double edge_length)
   {
        make_cone(CGAL::to_double(p0.x()), CGAL::to_double(p0.y()),CGAL::to_double(p0.z()),
                  CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),CGAL::to_double(p1.z()),
                  r0,r1,edge_length);
   }

  /**
   * \struct 
   * @brief Implicit function for a sphere with static center and radius.
   * 
   */
  struct sphere_wrapper{      
     public:
        static double radius;
        static double x0;
        static double y0;
        static double z0;
        static double function(double x, double y, double z){return  (x-x0)*(x -x0) +  (y-y0)*(y -y0) + (z-z0)*(z -z0) -radius*radius;}
        
   };
   
  // DocString: make_sphere
  /**
   * @brief Creates a sphere surface mesh.
   * 
   * Uses an implicit function to construct a structure of vertices and facets connecting vertices in the shape of 
   * a sphere. 
   *
   * @param x0 the x-coordinate of the center
   * @param y0 the y-coordinate of the center
   * @param z0 the z-coordinate of the center
   * @param r0 the radius  of  the sphere
   * @param mesh_resolution ratio between the sphere radius divided by the maximum edge length of the resulting surface mesh. 
   */
   template<typename Domain>
   void make_sphere( double x0, double y0, double  z0,double r0, double edge_length, double error_bound) 
   {    
      if( 2*CGAL_PI*r0<3*edge_length ) 
         throw InvalidArgumentError("Select smaller edge length."); 
  
      sphere_wrapper sphere;
      sphere.radius = r0;
      sphere.x0=x0;
      sphere.y0=y0;
      sphere.z0=z0;
      
     auto bounding = Kernel::Sphere_3(Point_3(x0,y0,z0), FT(r0*r0));   
     std::shared_ptr<Domain>   domain( new Domain(sphere.function, bounding,  error_bound));
     domain->create_surface_mesh(edge_length);
     auto surface = domain->template get_boundary<Surface>(1);
     this->mesh = surface->get_mesh();
     
     
   }
   
  // DocString: make_sphere   
  /**
   *  @brief Constructs a sphere surface.
   *  @param p0 SVMTK 3D point of the sphere center.
   *  @param r0 the radius of the sphere
   *  @param edge_length the upper edge length of the constructed triangulated surface 
   *  @overload 
   */
   template<typename Domain>
   void make_sphere(Point_3 p0, double r0, double edge_length, double error_bound)
   {
        make_sphere<Domain>(CGAL::to_double(p0.x()), CGAL::to_double(p0.y()),CGAL::to_double(p0.z()), 
                    r0, edge_length, error_bound);
   }
   
   //    
  /**
   * @brief Creates a circle in a plane, with the plane defined by a point and vector.
   *
   * @param point CGAL Point_3 object.  
   * @param vector CGAL Vector_3 object. 
   * @param radius of the circle.
   * @param edge_length the target edge_length for the resulting circle.
   */
   void make_circle_in_plane(Point_3 point, Vector_3 vector, double radius, double edge_length) 
   {
      face_vector fv1,fv3;
      Index v0 = mesh.add_vertex(point);
      Index vb, vt;
      Plane_3 plane(point,vector);  
      Vector_3 t1 = plane.base1();
      Vector_3 t2 = plane.base2();
      t1=t1/std::sqrt(t1.squared_length()); 
      t2=t2/std::sqrt(t2.squared_length()); 
      int number_of_segments = static_cast<int>(CGAL_PI*2*radius/edge_length); 
           
      if( number_of_segments<3 ) 
         throw InvalidArgumentError("Select smaller edge length.");      
         
      double c  = 360.0/(double)number_of_segments;
      for(int i = 0; i < number_of_segments; ++i) 
      {              
          Point_3 pb= point + t1*radius*std::cos(c*i*CGAL_PI/180) + t2*radius*std::sin(c*i*CGAL_PI/180); 
          vb = mesh.add_vertex(pb);
          if( i!=0 ) 
             fv1.push_back(mesh.add_face(v0,vb,Index(vb-1)));

      }
      fv1.push_back(mesh.add_face(Index(0), Index(1), vb));
      CGAL::Polygon_mesh_processing::isotropic_remeshing(fv1,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));
   }    
   
   // DocString: make_circle_in_plane  
   /**
   *  @brief Constructs a circle in the a given surface plane.
   *  @param px x-coordinate of point 
   *  @param py y-coordinate of point 
   *  @param pz z-coordinate of point
   *  @param vx x-coordinate of point 
   *  @param vy y-coordinate of point  
   *  @param vz z-coordinate of point             
   *  @param radius radius of circle 
   *  @param edge_length the upper edge length of the constructed triangulated surface 
   *  @overload
   */ 
   void make_circle_in_plane(double px, double py, double pz, double vx, double vy, double vz, double radius, double edge_length)
   {
        make_circle_in_plane(Point_3(px,py,pz),Vector_3(vx,vy,vz), radius, edge_length);
   } 
   
  // DocString: split_edges 
  /**  
   * @brief splits edges to a target edge length
   * 
   * CGAL function for splitting long edges.
   * @see [split_long_edges](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html)
   * @param target_edge_length the maximum edge legnth that the longer edges can split into.
   */
   void split_edges(double target_edge_length) 
   {
      assert_non_empty_mesh();
      CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
   }

  // DocString: save 
  /** 
   * @brief Saves the surface mesh to file.
   *
   * Valid file formats: off and stl. 
   * @param outpath string path to save file.
   *        @extensions : stl and off.
   */
   void save(const std::string outpath)
   {    
      assert_non_empty_mesh();  
      set_outward_face_orientation();  
      std::string extension = outpath.substr(outpath.find_last_of(".")+1);
      std::ofstream out(outpath);
      if( extension=="off" )
      {
         out << mesh;
         out.close();
      } 
      else if( extension=="stl" )
         CGAL::IO::write_STL(out,mesh);
   }
   
  // DocString: convex_hull
  /**
   * @brief Computes the convex hull of the points.
   * @see [CGAL::convex_hull_3](https://doc.cgal.org/latest/Convex_hull_3/index.html)
   * @returns a SVMTK Surface object.  
   */
   std::shared_ptr<Surface> convex_hull()
   {
      assert_non_empty_mesh(); 
      Polyhedron_3 poly;
      auto point_vector = get_points(); 
      CGAL::convex_hull_3(point_vector.begin(), point_vector.end(), poly);
      auto result = std::make_shared< Surface >(Surface(poly));
      return result;
   }

  /** 
   * @brief  Returns a vector that contains the point and the corresponding computed vertex normal.
   * @returns a vector that contains the point and the corresponding computed vertex normal.
   */
   std::vector<std::pair<Point_3,Vector_3>> get_points_with_normal()  
   {
      assert_non_empty_mesh();  
      std::vector<std::pair<Point_3, Vector_3>> result; 
      for(vertex_descriptor vit : mesh.vertices())
      { 
          Vector_3 normal = CGAL::Polygon_mesh_processing::compute_vertex_normal(vit,mesh);
          Point_3 point = mesh.point(vit); 
          result.push_back(std::make_pair(point,normal));
      }
      return result;
   } 
   
   struct find_largest_displacement{
         bool operator()(const std::pair<vertex_descriptor,Vector_3>&  a,  const std::pair<vertex_descriptor,Vector_3> & b) const 
                       { return( a.second.squared_length() > b.second.squared_length() );}
   };    

   /**
    * @brief Computes and returns the SVMTK Domain mesh resolution parameter of the surface.
    * @returns mesh resoltuion of the surface, ratio between minimum bounding radius and edge length. 
    */    
    double get_mesh_resolution()
    {
        double ael = average_edge_length();
        double r =  get_bounding_radius();
        return r/ael;
    }       
    
   /**
    * @brief Computes and returns the bounding radius, i.e. radius of a sphere that encloses the surface.
    * @returns bounding radius of the surface 
    */     
    double get_bounding_radius()
    {   
        auto bbox_3 = CGAL::Polygon_mesh_processing::bbox(mesh);
        
        double dx =  (bbox_3.max(0) - bbox_3.min(0))/2;
        double dy =  (bbox_3.max(1) - bbox_3.min(1))/2;        
        double dz =  (bbox_3.max(2) - bbox_3.min(2))/2; 
              
        return CGAL::sqrt(dx*dx+dy*dy+dz*dz); 
    }
    
   /**
    * @brief Returns AABB tree of the surface mesh.
    * @returns  AABB tree 
    */      
    AABB_Tree get_AABB_tree()
    {
       AABB_Tree tree(faces(mesh).first, faces(mesh).second, mesh);
       return tree;
    } 
    
   // DocString: get_collision_distance
   /**
     * @brief Computes the interesect distance in the facet normal direction. 
     * @returns vector of pairs, consistent of triangle points and distance  
     */       
    std::vector<std::pair<Triangle_3,double>> get_collision_distance(Surface other)
    {  
       std::vector<std::pair<Triangle_3,double>> Tridepth;
       double distance;       
     
       struct Skip
       {
           face_descriptor fd;
           Skip(const face_descriptor fd) : fd(fd) {}
           bool operator()(const face_descriptor& t) const
           { 
              return(t == fd);
           }
       };
       
       AABB_Tree tree = get_AABB_tree();
       AABB_Tree branch = other.get_AABB_tree();
      
       for(face_descriptor fd : faces(mesh))
       {
          halfedge_descriptor hd = halfedge(fd,mesh);

          Point_3 p1 = mesh.point(source(hd,mesh));
          Point_3 p2 = mesh.point(target(hd,mesh));
          Point_3 p3 = mesh.point(target(next(hd,mesh),mesh));

          Point_3 facetpoint = CGAL::centroid(p1,p2,p3);
          Vector_3 normal = CGAL::Polygon_mesh_processing::compute_face_normal(fd,mesh);
                    
          Ray ray(facetpoint, -normal);
          Skip skip(fd);
          Ray_intersection intersection = tree.first_intersection(ray, skip);
          distance = 0;
          
          if(intersection)
          {
             if(std::get_if<Point_3>(&(intersection->first)))
             {
                const Point_3* intersectpoint =  std::get_if<Point_3>(&(intersection->first));  
                if( branch.squared_distance(*intersectpoint) < FT(1.0e-18) )  
                    distance = CGAL::sqrt(CGAL::squared_distance(facetpoint, *intersectpoint));  
                else 
                    distance = -CGAL::sqrt(CGAL::squared_distance(facetpoint, *intersectpoint));  
             }
         }
         Triangle_3 tri(p1,p2,p3);
         Tridepth.push_back(std::make_pair(tri,distance));    
       }   
       return Tridepth;
    }
    
   // DOCSTRING get_collision_spheres 
   /**
    * 
    * @brief Computes the collision sphere radius for each facet in the normal direction.
    *  
    * @returns  
    */       
    std::vector<std::pair<Triangle_3,double>> get_collision_spheres()
    {  
        std::vector<std::pair<Triangle_3,double>> Tridepth;    
        double beta,min_beta;
        
        double proximity = std::min(1000, num_vertices() );  
        
        Vertex_point_pmap vppmap = get(CGAL::vertex_point,mesh); 
     
        Tree tree(vertices(mesh).begin(),
                 vertices(mesh).end(),
                 Splitter(),
                 Traits(vppmap));

        Distance tr_dist(vppmap);
        for(face_descriptor fd : faces(mesh)) 
        {
          min_beta = 100.0; // Large value
          halfedge_descriptor hd = halfedge(fd,mesh);

          Point_3 p1 = mesh.point(source(hd,mesh));
          Point_3 p2 = mesh.point(target(hd,mesh));
          Point_3 p3 = mesh.point(target(next(hd,mesh),mesh));
          
          Point_3 facetpoint = CGAL::centroid(p1,p2,p3);
          
          Vector_3 normal = CGAL::Polygon_mesh_processing::compute_face_normal(fd,mesh);

          K_neighbor_search search(tree,facetpoint,  proximity, 0, true, tr_dist); 
          
          for( auto pit = search.begin() ; pit != search.end() ; pit++) 
          {
             Vector_3 b(mesh.point(pit->first), facetpoint);
             // Pythagorean theorem : c^2 = b^2 + a^2 -2 a*b , -> c = a = -n beta
             if( CGAL::abs(-normal*b) < FT(0.00001) )  
                 beta =  abs(b.squared_length()/(0.00001*2));                    
             else  
                 beta =  abs(b.squared_length()/(-normal*b*2));  
              
             if( beta<min_beta ) 
                 min_beta=beta;
          }  
          Triangle_3 tri(p1,p2,p3);
          Tridepth.push_back(std::make_pair(tri,min_beta));
       }
       return Tridepth;
    }  

   /**
    * @brief Computes 
    * @param none 
    * @return the avearge triangle area. 
    */
    double average_triangle_area()
    {
        assert_non_empty_mesh();
        double sum_area;
        Point_3 p1,p2,p3;

        for(face_descriptor fd : faces(mesh)) 
        {    
             halfedge_descriptor hd = halfedge(fd,mesh);
             p1 = mesh.point(source(hd,mesh));
             p2 = mesh.point(target(hd,mesh));
             p3 = mesh.point(target(next(hd,mesh),mesh));
             sum_area += CGAL::to_double( CGAL::sqrt(Triangle_3(p1,p2,p3).squared_area()) );
        }
        return sum_area/num_faces(); 
    }
    
    
    /**
     * @brief Finds "degenerate" faces
     * @param none
     * @return vector of faces 
     */
    face_vector  get_degenerated_faces()
    {
         assert_non_empty_mesh();
         Point_3 p1,p2,p3;
         face_vector output;

         FT avg_edge(average_edge_length());
         FT avg_area(average_triangle_area());
         
         for(face_descriptor fd : faces(mesh)) 
         {
             halfedge_descriptor hd = halfedge(fd,mesh);
             p1 = mesh.point(source(hd,mesh));
             p2 = mesh.point(target(hd,mesh));
             p3 = mesh.point(target(next(hd,mesh),mesh));
             
             if (  Triangle_3(p1,p2,p3).squared_area() < 0.04*avg_area*avg_area)
             {   
                output.push_back(fd);
                continue;
             }
           
             if ( CGAL::squared_distance(p1,p2) < 0.04*avg_edge*avg_edge or
                  CGAL::squared_distance(p2,p3) < 0.04*avg_edge*avg_edge or 
                  CGAL::squared_distance(p3,p1) < 0.04*avg_edge*avg_edge) 
             {
                   output.push_back(fd);
             }
         }
         return output;
    }
    
    // 
    /**
     * @brief Finds and removes "degenerate" faces
     * @param none 
     * @return
     */
    void remove_degenerate_faces()
    {
        auto removal = get_degenerated_faces();
        remove_faces(removal);         
    }
   
    Vector_3 average_normal_vertex( vertex_vector&  vertices)
    {
        Vector_3 n1 ;
        Vector_3 delta=CGAL::NULL_VECTOR;
        for( auto vit : vertices) 
        {
             delta +=  CGAL::Polygon_mesh_processing::compute_vertex_normal(vit,mesh);   
        }    
        n1= delta/CGAL::sqrt(delta.squared_length());
        return n1;
    }
    
   /**
    * @brief Sets the proximity ratio. 
    * @note The proximity ratio is used to determine close vertice. 
    *       If the distance between two vertices exceeds the proximity ratio times edge length,
    *       then the vertices are considered close.   
    *       Values should be between 1 and 2. 
    * @param proximity_ratio 
    */
    void set_proximity_ratio(double proximity_ratio ) 
    {
         _prratio = proximity_ratio; 
    }

   /**
    * @brief Sets the smoothing reduction. 
    * @note When the displacment is lower than the lower displacment bound, 
    *       then the smoothing reduction is used to decrease the smoothing.  
    */
    void set_smoothing_reduction(double smooth_reduction) 
    {
        _smreduc = smooth_reduction;
    }
            
   /**
    * @brief Sets the displacement ratio. If the maximum displacement is lower 
    * @note  The displacment ratio times average egde length determines the lower displacement bound 
    * 
    */ 
    void set_displacment_ratio(double displacment_ratio) 
    {
          _stratio = displacment_ratio;
    }

   /**
    * @brief  
    * 
    * 
    */
    double get_smoothing_reduction_factor() 
    {
           return _smreduc;
    }

   /**
    * @brief 
    * 
    * 
    */       
    double get_lower_displacment_bound() 
    {
       return average_edge_length()*average_edge_length()*_stratio;          
    }  
  
    
   protected:
    Mesh mesh;
    
    double  _stratio; 
    double  _prratio; 
    double  _smreduc;     
};

/**
 * @brief Inline declaration of static variables in header files requrie c++ version equal or greate than 7 
 *       with the compiler option c++1z/c++17.
 * @hideinitializer 
*/
inline double Surface::sphere_wrapper::radius = 0; // inline declaration of static variable radius for sphere wrapper.
inline double Surface::sphere_wrapper::x0 = 0;     // inline declaration of static variable radius for sphere wrapper.
inline double Surface::sphere_wrapper::y0 = 0;     // inline declaration of static variable radius for sphere wrapper.
inline double Surface::sphere_wrapper::z0 = 0;     // inline declaration of static variable radius for sphere wrapper.

/** 
 * @brief Constructs a convex hull from a vector of points
 *
 * Uses CGAL convex hull algorithm on a vector of points. 
 *
 * @param point_vector a vector of template points, default option is surface points.
 * @returns SVMTK Surface object.
 */
template< typename Surface, typename Point_3 = typename Surface::Point_3>
std::shared_ptr<Surface> convex_hull(std::vector<Point_3 >& point_vector)
{
  typename Surface::Polyhedron_3 polyhedron;
  CGAL::convex_hull_3(point_vector.begin(), point_vector.end(), polyhedron);
  return std::make_shared<Surface>(Surface(polyhedron));
}

// DocString: separate_surface_overlapp
/** 
 * @brief Separates two overlapping surfaces by contraction of surfaces boundary.
 *
 * Separates two overlapping surfaces iteratively by moving the vertices in the negative normal 
 * direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each 
 * vertex. This continues until the closest vertex is adjacent.
 * The centeroids of the surfaces should be sufficiently apart, compared to volume.   
 *
 * @param surf1 first SVMTK Surface object.  
 * @param surf2 second SVMTK Surface object. 
 * @param edge_movement the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param smoothing the input of the laplacian smoothing that is applied after each iteration of vertex movement. 
 * @param max_iter maximum number of iteration.
 * @returns true if complete.
 */
template<typename Surface>
bool separate_surface_overlapp(Surface& surf1, Surface& surf2, double edge_movement, double smoothing, int max_iter)
{
  typedef typename Surface::vertex_vector_map vertex_vector_map;
  typedef typename Surface::vertex_vector vertex_vector;

  vertex_vector s1vertices  = surf1.get_vertices_inside(surf2);
  vertex_vector s2vertices  = surf2.get_vertices_inside(surf1);

  double s1ael = surf1.get_lower_displacment_bound(); //average_edge_length();
  double s2ael = surf2.get_lower_displacment_bound(); //average_edge_length();

  double srf1 = surf1.get_smoothing_reduction_factor();
  double srf2 = surf2.get_smoothing_reduction_factor();
   //get_lower_displacment_bound() 

  double smth1 = smoothing;
  double smth2 = smoothing;   
  
  int iter =0;
  while( !s1vertices.empty() or !s2vertices.empty() ) 
  {      
        vertex_vector_map vertex_displacement_1 = surf1.get_vertex_displacement(s1vertices, -abs(edge_movement),smth1);     
        vertex_vector_map vertex_displacement_2 = surf2.get_vertex_displacement(s2vertices, -abs(edge_movement),smth2);         

        surf1.set_adjacent_vertices(vertex_displacement_1);
        surf2.set_adjacent_vertices(vertex_displacement_2);    
      
        auto max_element_1 = std::max_element( vertex_displacement_1.begin(), vertex_displacement_1.end(), typename Surface::find_largest_displacement());    
        auto max_element_2 = std::max_element( vertex_displacement_2.begin(), vertex_displacement_2.end(), typename Surface::find_largest_displacement());    

        if(  max_element_1->second.squared_length() < s1ael ) 
             smth1 = srf1*smth1; 
        if ( max_element_2->second.squared_length() < s2ael ) 
             smth2 = srf2*smth2;
             
        if( surf1.num_self_intersections()>0  )
        {
            surf1.repair_self_intersections();
            s1vertices  = surf1.get_vertices_inside(surf2);
            smth1 = smth1;       
        }
        else 
            s1vertices  = surf1.get_vertices_inside(surf2,s1vertices);
        
        if( surf2.num_self_intersections()>0 )
        {
            surf2.repair_self_intersections();
            s2vertices  = surf2.get_vertices_inside(surf1);
            smth2 = smth2;     
        }
        else 
            s2vertices  = surf2.get_vertices_inside(surf1,s2vertices);
   
       if( s1vertices.empty() )  
        {
            s1vertices = surf1.get_vertices_inside(surf2);
            smth1 = smth1;    
        }             

        if( s2vertices.empty() )  
        {
            s2vertices = surf2.get_vertices_inside(surf1);
            smth2 = smth2;    
        }
  
       if( iter++>max_iter )
       {
         std::cout << "Failed to converge in "<< max_iter <<" steps, terminating" << std::endl;
         return false;
       }
       
  }
  return true;
}

// DocString: separate_surface_overlapp
/** 
 * @brief Separates two overlapping surfaces outside a thrid surface.
 *
 * Separates the points of two surfaces that are outside a third surface iteratively by moving the vertices in the negative normal 
 * direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each 
 * vertex. This continues until the closest vertex is adjacent.
 *
 * @param surf1 first SVMTK Surface object.  
 * @param surf2 second SVMTK Surface object. 
 * @param other third SVMTK Surface object, the algorithm is not applied to vertices inside this surface.
 * @param edge_movement the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param smoothing the input of the laplacian smoothing that is applied after each iteration of vertex movement. 
 * @param max_iter maximum number of iteration.
 * @returns true if complete.
 */
template<typename Surface>
bool separate_surface_overlapp(Surface& surf1, Surface& surf2, Surface& other, double edge_movement, double smoothing, int max_iter)
{
   typedef typename Surface::vertex_vector_map vertex_vector_map;
   typedef typename Surface::vertex_vector vertex_vector;

   vertex_vector s1vertices  = surf1.get_vertices_inside(surf2);  
   vertex_vector s2vertices  = surf2.get_vertices_inside(surf1);

   double s1ael = surf1.get_lower_displacment_bound(); //average_edge_length();
   double s2ael = surf2.get_lower_displacment_bound(); //average_edge_length();

   double srf1 = surf1.get_smoothing_reduction_factor();
   double srf2 = surf2.get_smoothing_reduction_factor();
  
   int iter =0;   
   double smth1 = smoothing;
   double smth2 = smoothing;   
    
   while( !s1vertices.empty() or !s2vertices.empty() ) 
   {
        s1vertices  = surf1.get_vertices_outside(other,s1vertices); 
        s2vertices  = surf2.get_vertices_outside(other,s2vertices);

        vertex_vector_map vertex_displacement_1 = surf1.get_vertex_displacement(s1vertices, -abs(edge_movement), smth1);     
        vertex_vector_map vertex_displacement_2 = surf2.get_vertex_displacement(s2vertices, -abs(edge_movement), smth2);         

        surf1.set_adjacent_vertices(vertex_displacement_1);
        surf2.set_adjacent_vertices(vertex_displacement_2);           
        
        auto max_element_1 = std::max_element( vertex_displacement_1.begin(), vertex_displacement_1.end(), typename Surface::find_largest_displacement());    
        auto max_element_2 = std::max_element( vertex_displacement_2.begin(), vertex_displacement_2.end(), typename Surface::find_largest_displacement());    

        if(  max_element_1->second.squared_length() < s1ael ) 
             smth1 = srf1*smth1; 
        if ( max_element_2->second.squared_length() < s2ael ) 
             smth2 = srf2*smth2;
             
                     
        if( surf1.num_self_intersections()>0 ) 
        {
            surf1.repair_self_intersections();
            s1vertices = surf1.get_vertices_inside(surf2);
            smth1 = smoothing;    
        }
        else 
            s1vertices = surf1.get_vertices_inside(surf2,s1vertices);
        
        if( surf2.num_self_intersections()>0 )
        {
            surf2.repair_self_intersections();
            s2vertices = surf2.get_vertices_inside(surf1);
            smth2 = smoothing;      
        }
        else 
            s2vertices = surf2.get_vertices_inside(surf1,s2vertices);

        if( s1vertices.empty() )  
        {
            s1vertices = surf1.get_vertices_inside(surf2);
            smth1 = smoothing;    
        }             

        if( s2vertices.empty() )  
        {
            s2vertices = surf2.get_vertices_inside(surf1);
            smth2 = smoothing;    
        }
    
        s1vertices  = surf1.get_vertices_outside(other,s1vertices); 
        s2vertices  = surf2.get_vertices_outside(other,s2vertices);

        if( iter++>max_iter)
        {
           std::cout << "Failed to converge in "<< max_iter <<" steps, terminating" << std::endl;
           return false;       
        }
   }  
   return true;
}

// DocString: separate_close_overlapp
/**
 * @brief Separates two close surfaces outside a thrid surface.
 *
 * Separates the points of two surfaces that are outside a third surface iteratively by moving the vertices in the negative normal 
 * direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each 
 * vertex.This continues until the closest vertex is adjacent.
 *
 * The implementation uses -edge_movement for backwards compatbility in the documentation.
 * 
 * @param surf1 first SVMTK Surface object.  
 * @param surf2 second SVMTK Surface object. 
 * @param other third SVMTK Surface object, the algorithm is not applied to vertices inside this surface.
 * @param edge_movement the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param smoothing the input of the laplacian smoothing that is applied after each iteration of vertex movement 
 * @param max_iter maximum number of iteration.
 * @returns true if complete. 
 * @note  for backwards compatiblity : abs(edge_movement) 
 */
template< typename Surface>
bool separate_close_surfaces(Surface& surf1, Surface& surf2, Surface& other, double edge_movement, double smoothing, int max_iter)
{
   typedef typename Surface::vertex_vector_map vertex_vector_map;
   typedef typename Surface::vertex_vector vertex_vector;

   vertex_vector s1vertices =  surf1.get_close_vertices(surf2);
   vertex_vector s2vertices =  surf2.get_close_vertices(surf1);

   
   double s1ael = surf1.get_lower_displacment_bound(); //average_edge_length();
   double s2ael = surf2.get_lower_displacment_bound(); //average_edge_length();

   double srf1 = surf1.get_smoothing_reduction_factor();
   double srf2 = surf2.get_smoothing_reduction_factor();
   
   
   int iter=0;
   double smth1 = smoothing;
   double smth2 = smoothing;
  
   s1vertices  = surf1.get_vertices_outside(other,s1vertices);
   s2vertices  = surf2.get_vertices_outside(other,s2vertices);      
  
   while( !s1vertices.empty() or !s2vertices.empty() ) 
   {
        vertex_vector_map vertex_displacement_1 = surf1.get_vertex_displacement(s1vertices, -abs(edge_movement), smth1);     
        vertex_vector_map vertex_displacement_2 = surf2.get_vertex_displacement(s2vertices, -abs(edge_movement), smth2);         
 
        auto max_element_1 = std::max_element( vertex_displacement_1.begin(), vertex_displacement_1.end(), typename Surface::find_largest_displacement());    
        auto max_element_2 = std::max_element( vertex_displacement_2.begin(), vertex_displacement_2.end(), typename Surface::find_largest_displacement());    
        
        surf1.set_adjacent_vertices(vertex_displacement_1);
        surf2.set_adjacent_vertices(vertex_displacement_2);           
        
        if(  max_element_1->second.squared_length() < s1ael ) 
             smth1 = srf1*smth1; 
        if ( max_element_2->second.squared_length() < s2ael ) 
             smth2 = srf2*smth2;
        
        if( surf1.num_self_intersections()>0 )
        {
            surf1.repair_self_intersections();
            s1vertices = surf1.get_close_vertices(surf2); 
            s1vertices  = surf1.get_vertices_outside(other,s1vertices);  
            smth1 = smoothing;    
        }
        else 
            surf1.get_close_vertices(surf2,s1vertices); 
   
        if( surf2.num_self_intersections()>0   )
        {
            surf2.repair_self_intersections();  
            s2vertices = surf2.get_close_vertices(surf1); 
            s2vertices = surf2.get_vertices_outside(other,s2vertices);    
            smth2 = smoothing;   
        }
        else 
            surf2.get_close_vertices(surf1,s2vertices); 


        if( s1vertices.empty() )  
        {
            s1vertices = surf1.get_close_vertices(surf2); 
            s1vertices = surf1.get_vertices_outside(other,s1vertices);   
            smth1 = smoothing;    
        }             

        if( s2vertices.empty() )  
        {
            s2vertices = surf2.get_close_vertices(surf1); 
            s2vertices = surf2.get_vertices_outside(other,s2vertices);   
            smth2 = smoothing;    
        }

        if( iter++>max_iter ) 
        {
           std::cout << "Failed to converge in "<< max_iter <<" steps, terminating" << std::endl;
           return false;                
        }
   }

   return true;
}

// DocString: separate_close_overlapp
/**
 * @brief Separates two close surfaces.
 * 
 * Separates two surfaces iteratively by moving the vertices in the negative normal 
 * direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each 
 * vertex. This continues until the closest vertex is adjacent. 
 *
 * The implementation uses -edge_movement for backwards compatbility.
 * 
 * @param surf1 first SVMTK Surface object.  
 * @param surf2 second SVMTK Surface object. 
 * @param edge_movement the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param smoothing the input of the laplacian smoothing that is applied after each iteration of vertex movement. 
 * @param max_iter maximum number of iteration.
 * @returns true if complete.
 */
template<typename Surface> 
bool separate_close_surfaces(Surface& surf1, Surface& surf2, double edge_movement, double smoothing, int max_iter)
{
   typedef typename Surface::vertex_vector_map vertex_vector_map;
   typedef typename Surface::vertex_vector vertex_vector;  
 
   vertex_vector s1vertices = surf1.get_close_vertices(surf2);   
   vertex_vector s2vertices = surf2.get_close_vertices(surf1);

   double s1ael = surf1.get_lower_displacment_bound(); //average_edge_length();
   double s2ael = surf2.get_lower_displacment_bound(); //average_edge_length();

   double srf1 = surf1.get_smoothing_reduction_factor();
   double srf2 = surf2.get_smoothing_reduction_factor();
   
   double smth1 = smoothing;
   double smth2 = smoothing;
   
   int iter=0;   
   while( !s1vertices.empty() or !s2vertices.empty() ) 
   {
        vertex_vector_map vertex_displacement_1 = surf1.get_vertex_displacement(s1vertices, -abs(edge_movement), smth1);     
        vertex_vector_map vertex_displacement_2 = surf2.get_vertex_displacement(s2vertices, -abs(edge_movement), smth2);         

        surf1.set_adjacent_vertices(vertex_displacement_1);
        surf2.set_adjacent_vertices(vertex_displacement_2);  
        
        auto max_element_1 = std::max_element( vertex_displacement_1.begin(), vertex_displacement_1.end(), typename Surface::find_largest_displacement());    
        auto max_element_2 = std::max_element( vertex_displacement_2.begin(), vertex_displacement_2.end(), typename Surface::find_largest_displacement());   

        if(  max_element_1->second.squared_length() < s1ael ) 
             smth1 = srf1*smth1; 
        if ( max_element_2->second.squared_length() < s2ael ) 
             smth2 = srf2*smth2;
                    
        if( surf1.num_self_intersections()>0 )
        {
            surf1.repair_self_intersections();
            s1vertices = surf1.get_close_vertices(surf2); 
            smth1 = smoothing;      
        }
        else 
            surf1.get_close_vertices(surf2,s1vertices); 
   
        if( surf2.num_self_intersections()>0 )
        {
            surf2.repair_self_intersections();  
            s2vertices = surf2.get_close_vertices(surf1); 
            smth2 = smoothing;      
        }
        else 
            surf2.get_close_vertices(surf1,s2vertices);   
        
        
        if( s1vertices.empty() )  
        {
            s1vertices = surf1.get_close_vertices(surf2); 
            smth1 = smoothing;    
        }    
            
        if( s2vertices.empty() )  
        {
            s2vertices = surf2.get_close_vertices(surf1); 
            smth2 = smoothing;       
        }

        if( iter++>max_iter ) 
        {
           std::cout << "Failed to converge in "<< max_iter <<" steps, terminating" << std::endl;
           return false;                
        }
        
   } 
   return true;
}


// DocString: union_partially_overlapping_surfaces
/**
 * @brief Takes surface union of surfaces that partially overlapp each other. 
 * 
 * Requires that the surfaces does overlap. 
 * Surface vertices inside the other surface are obtained. Then adjacent vertices
 * within a threshold angle of the vertex normal are added to the vertex vector.
 * The points corresponding to the vertex vector are iteratively moved in the 
 * vertex normal direction until there is suffient overlapp between surfaces.             
 *
 * @tparam SVMTK Surface class.
 * @param SVMTK Surface object.  
 * @param SVMTK Surface object.
 * @param angle_in_degree 
 * @param adjustment the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param smoothing the number of taubin iterations after each iteration of vertex movement.
 * @param max_iter maximum number of iteration. 
 * @returns a SVMTK Surface object.
 */
template<typename Surface> 
std::shared_ptr<Surface> union_partially_overlapping_surfaces( Surface& surf1, Surface& surf2, double angle_in_degree, double adjustment, double smoothing, int max_iter )
{
     typedef typename Surface::vertex_vector_map vertex_vector_map;
     typedef typename Surface::vertex_vector vertex_vector;
     typedef typename Surface::Vector_3 Vector_3; 

     vertex_vector s1vertices,s2vertices;
     
     s1vertices  = surf1.get_vertices_inside(surf2); 
     s2vertices  = surf2.get_vertices_inside(surf1);

     Vector_3 disp1 = adjustment*surf1.average_normal_vertex(s1vertices);
     Vector_3 disp2 = adjustment*surf2.average_normal_vertex(s2vertices);      

     surf1.add_adjacent_vertices(s1vertices,angle_in_degree);
     surf2.add_adjacent_vertices(s2vertices,angle_in_degree);    

     double smth1 = smoothing;
     double smth2 = smoothing;
    
     double avg_edge_length = 0.5*(surf1.average_edge_length() + surf2.average_edge_length());
    
     double s1ael = surf1.get_lower_displacment_bound(); //average_edge_length();
     double s2ael = surf2.get_lower_displacment_bound(); //average_edge_length();

     double srf1 = surf1.get_smoothing_reduction_factor();
     double srf2 = surf2.get_smoothing_reduction_factor();
   
     int iter =0;
     while( !s1vertices.empty() or !s2vertices.empty() )   
     {     
   
           vertex_vector_map vertex_displacement_1 = surf1.get_constant_vertex_displacement(s1vertices, disp1 , smth1);
           vertex_vector_map vertex_displacement_2 = surf2.get_constant_vertex_displacement(s2vertices, disp2 , smth2);                
            
           auto max_element_1 = std::max_element( vertex_displacement_1.begin(), vertex_displacement_1.end(), typename Surface::find_largest_displacement());    
           auto max_element_2 = std::max_element( vertex_displacement_2.begin(), vertex_displacement_2.end(), typename Surface::find_largest_displacement());    
                  

           surf1.set_adjacent_vertices(vertex_displacement_1); 
           surf2.set_adjacent_vertices(vertex_displacement_2);            
                  
           if(  max_element_1->second.squared_length() < s1ael ) 
               smth1 = srf1*smth1; 
           if ( max_element_2->second.squared_length() < s2ael ) 
               smth2 = srf2*smth2;
          
           if( s1vertices.empty()) 
           {
               s1vertices  = surf1.get_vertices_inside(surf2); 
               surf1.add_adjacent_vertices(s1vertices,angle_in_degree);
               s1vertices = surf1.get_vertices_outside(surf2,s1vertices);   
               smth1 = smoothing;          
           }
           
           if( s2vertices.empty() ) 
           {
               s2vertices  = surf2.get_vertices_inside(surf1);
               surf2.add_adjacent_vertices(s2vertices,angle_in_degree);
               smth2 = smoothing;         
           } 

    
          if( surf1.num_self_intersections()>0 )
          {
            surf1.repair_self_intersections();
            s1vertices  = surf1.get_vertices_inside(surf2); 
            surf1.add_adjacent_vertices(s1vertices,angle_in_degree);
            smth1 = smoothing;      
          }
     
          if( surf2.num_self_intersections()>0 )
          {
             surf2.repair_self_intersections();  
             s1vertices  = surf1.get_vertices_inside(surf2); 
             surf1.add_adjacent_vertices(s1vertices,angle_in_degree);
             smth2 = smoothing;      
          }

           s1vertices = surf1.get_vertices_outside(surf2,s1vertices); 
           s2vertices = surf2.get_vertices_outside(surf1,s2vertices);

           if( iter++>max_iter )
               break;
     
     }

     surf1.repair_self_intersections();
     surf2.repair_self_intersections();     
     
     std::shared_ptr<Surface> result(new Surface()); 
     try
     {
       CGAL::Polygon_mesh_processing::corefine_and_compute_union(surf1.get_mesh(), surf2.get_mesh(), result->get_mesh());     
     }
     catch (const std::exception &exc)
     {
       std::string output = "CGAL precondition error\n"+ surf1.CGAL_precondition_evaluation(surf2);
       throw PreconditionError(output.c_str());
     }  
       
     result->smooth_taubin(20);
     result->isotropic_remeshing(0.5*avg_edge_length, 5, false);


     return result;
}

template<typename Surface> 
bool enclose( Surface& surf1, Surface& surf2, double adjustment, double smoothing, int max_iter )
{
     auto ael = surf1.average_edge_length();
     auto result_1 = surf1.enclose(surf2 , abs(adjustment), smoothing, max_iter);
     surf1.isotropic_remeshing(ael,5,false);    
     surf1.smooth_taubin(5);
     auto result_2 = surf1.separate(surf2, abs(adjustment), smoothing, max_iter);
     surf1.isotropic_remeshing(ael,5,false);
     surf1.smooth_taubin(5);     
     return (result_1.first and result_2.first);
}

template<typename Surface> 
bool expose( Surface& surf1, Surface& surf2, double adjustment, double smoothing, int max_iter )
{
      auto ael = surf1.average_edge_length();
      auto result_1 =  surf1.expose(surf2  , -abs(adjustment),smoothing, max_iter);
      surf1.isotropic_remeshing(ael,5,false);
      surf1.smooth_taubin(5);     
      auto result_2 =  surf1.separate(surf2, -abs(adjustment),smoothing, max_iter);
      surf1.isotropic_remeshing(ael,5,false);
      surf1.smooth_taubin(5);      
      return (result_1.first and result_2.first);
}

template<typename Surface> 
bool embed( Surface& surf1, Surface& surf2, double adjustment, double smoothing, int max_iter )
{  
     auto ael = surf1.average_edge_length();
     auto result_1 = surf1.embed(surf2   , -abs(adjustment), smoothing, max_iter);
     surf1.isotropic_remeshing(ael,5,false);
     surf1.smooth_taubin(5);
     auto result_2 = surf1.separate(surf2, -abs(adjustment), smoothing, max_iter);
     surf1.isotropic_remeshing(ael,5,false);
     surf1.smooth_taubin(5);     
     return (result_1.first and result_2.first);
}





#endif
