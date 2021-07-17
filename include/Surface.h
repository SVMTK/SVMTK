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

#define PI 3.14159265

/* --Includes -- */
#include "surface_mesher.h"
#include "Errors.h"

/* -- boost-- */
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/system/error_code.hpp>

/* -- CGAL 2D and 3D Linear Geometry Kernel  -- */
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>

/* -- CGAL 3D Convex Hulls  -- */
#include <CGAL/convex_hull_3.h>

/* -- CGAL Principal Component Analysis -- */
#include <CGAL/centroid.h>

/* -- CGAL dD Spatial Searching  -- */
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

/* -- CGAL Triangulated Surface Mesh Shortest Paths-- */
#include <CGAL/Surface_mesh_shortest_path.h>

/* -- CGAL Triangulated Surface Mesh Skeletonization -- */
#include <CGAL/extract_mean_curvature_flow_skeleton.h>

/* -- CGAL IO -- */
#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/STL_writer.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>

/* -- CGAL Polygon Mesh Processing -- */
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
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

/* -- CGAL Surface mesh simplification v5.0.4 -- */
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>

/* -- CGAL AABB -- */
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

/* -- CGAL Triangulated Surface Mesh Segmentation -- */
#include <CGAL/mesh_segmentation.h>


//#include <CGAL/internal/Mesh_3/helpers.h>
//#include <CGAL/boost/graph/split_graph_into_polylines.h>
//#include <sys/stat.h>

// FIXME INLUCDE @throws

/**
 * \class 
 * Cost function used to collapse edges with 
 * zero cost of removal, like edges in a 
 * straight line 
 */
template<class TM_>
class Cost_stop_predicate
{
  public:
    typedef TM_ TM;

    Cost_stop_predicate(double threshold) : thres(threshold) {}
    /**
     * TODO
     */
    template< typename F, typename Profile>
    bool operator()( F const & aCurrentCost, Profile const & profile, std::size_t aic ,std::size_t acc) const
    {
      return static_cast<double>(aCurrentCost) > thres;
    }
  private :
    double thres;
};


/**
 * @brief Loads triangulated surfaces with exentsion off or stl.
 *
 * Loads surface mesh as a polygon soup, and 
 * reassembles it according to CGAL structure
 * of faces and vertices.
 * 
 * @param[in] file string of the fileaneme 
 * @param[out] mesh triangulated surface mesh.
 * @return true if loaded  
 */
template< typename Mesh >
bool load_surface(const std::string file, Mesh& mesh)
{
  typedef typename Mesh::Point Point_3;

  std::vector<Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;

  std::ifstream input(file);

  std::string extension = file.substr(file.find_last_of(".") + 1);
  if (!input) //FIXME throw 
  {
    throw InvalidArgumentError("Cannot open file");
    return false;
  }

  if (extension == "off")
  {
    if (!CGAL::read_OFF(input, points, polygons))
    {
      throw InvalidArgumentError("Error parsing the OFF file."); 
      return false;
    }
  }
  else if (extension == "stl")
  {
    if (!CGAL::read_STL(input, points, polygons))
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
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,mesh);
  CGAL::Polygon_mesh_processing::orient(mesh);

  if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
  {
    CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
  }
  return true;
}

/** 
 * @brief Constructs a convex hull from a vector of points
 *
 * Uses CGAL convex hull algorithm on a vector of points. TODO
 *
 * @param point_vector a vector of template points, default option is surface points 
 * @return a object of the Surface class
 */
template< typename Surface, typename Point_3 = typename Surface::Point_3>
std::shared_ptr< Surface > convex_hull(std::vector<Point_3 >& point_vector)
{
  typename Surface::Polyhedron poly;
  CGAL::convex_hull_3(point_vector.begin(), point_vector.end(), poly);
  auto result = std::make_shared< Surface >(Surface(poly));
  return result;
}

/** 

 * @brief Separates two overlapping surfaces by contraction of surfaces boundary.
 *
 * Separates two overlapping surfaces iteratively by moving the vertices in the negative normal 
 * direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each 
 * vertex. This continues until the closest vertex is adjacent.
 * The centeroids of the surfaces should be sufficiently apart, compared to volume.   
 *
 * @param surf1 first SVMTK surface class object  
 * @param surf2 second SVMTK surface class object 
 * @param double edge_movement the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param double smoothing the input of the laplacian smoothing that is applied after each iteration of vertex movement 
 * @param max_iter maximum number of iteration.
 * @return Surface class union of the modified 
 */
template< typename Surface>
bool separate_surface_overlapp(Surface& surf1 , Surface& surf2, double edge_movement=-0.25, double smoothing=0.3, int max_iter=400)
{
  typedef typename Surface::vertex_vector vertex_vector;
  typedef typename Surface::vertex_descriptor vertex_descriptor;

  double surf1_ael = surf1.average_edge_length(); // surf1.vertices_inside new ??
  double surf2_ael = surf2.average_edge_length();
  
  std::map< vertex_descriptor, double> map1,map2;
  vertex_vector surface1_vertices,surface2_vertices;
  surface1_vertices  = surf1.get_vertices_inside(surf2);
  surface2_vertices  = surf2.get_vertices_inside(surf1);

  int iter =0;

  while (!surface1_vertices.empty() or !surface2_vertices.empty())
  {
  
        surf1.get_adjacent_vertices(surface1_vertices);
        surf2.get_adjacent_vertices(surface2_vertices);
    
        map1 = surf1.get_shortest_edge_map(surface1_vertices,edge_movement);
        map2 = surf2.get_shortest_edge_map(surface2_vertices,edge_movement);

        surf1.adjust_vertices_in_region(map1.begin(),map1.end());
        surf2.adjust_vertices_in_region(map2.begin(),map2.end());
  
        surf1.smooth_laplacian_region(surface1_vertices.begin(), surface1_vertices.end(), smoothing);
        surf2.smooth_laplacian_region(surface2_vertices.begin(), surface2_vertices.end(), smoothing);

        if (surf1.num_self_intersections()>0)
        {
            surf1.collapse_edges(surf1_ael);
            surface1_vertices  = surf1.get_vertices_inside(surf2);
        
        }
        else 
            surface1_vertices  = surf1.get_vertices_inside(surf2,surface1_vertices);
        
        if (surf2.num_self_intersections()>0)
        {
            surf2.collapse_edges(surf2_ael);
            //surf2.isotropic_remeshing(surf2_ael, 1, false);
            surface2_vertices  = surf2.get_vertices_inside(surf1);
        }
        else 
            surface2_vertices  = surf2.get_vertices_inside(surf1,surface2_vertices);
    
    //surface1_vertices  = surf1.get_vertices_inside(surf2);
    //surface2_vertices  = surf2.get_vertices_inside(surf1);


    if ( surf1.num_self_intersections()+surf2.num_self_intersections() >0 )
    {
          std::cout << "Detected "<< surf1.num_self_intersections()+surf2.num_self_intersections() <<" self-intersections" << std::endl;
          std::cout << "Recommended to use istropic_remeshing before continuation"  << std::endl; // soft error handling ?? 
          return false;
    }
    if ( iter++>max_iter)
    {
      std::cout << "Failed to converge in "<< max_iter <<" steps, terminating" << std::endl;
      return false;
    }
  }
  return true;
}

/** 
 * @brief Separates two overlapping surfaces outside a thrid surface.
 *
 * Separates the points of two surfaces that are outside a third surface iteratively by moving the vertices in the negative normal 
 * direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each 
 * vertex. This continues until the closest vertex is adjacent.
 *
 * @param surf1 first SVMTK surface class object  
 * @param surf2 second SVMTK surface class object 
 * @param other third SVMTK surface class obejct, the algorithm is not applied to vertices inside this surface.
 * @param double edge_movement the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param double smoothing the input of the laplacian smoothing that is applied after each iteration of vertex movement 
 * @param max_iter maximum number of iteration.
 * @return Surface class union of the modified 
 */
template< typename Surface> 
bool separate_surface_overlapp(Surface& surf1 , Surface& surf2 , Surface& other, double edge_movement=-0.25, double smoothing=0.3, int max_iter=400)
{
    typedef typename Surface::vertex_vector vertex_vector;
    typedef typename Surface::vertex_descriptor vertex_descriptor;


    double surf1_ael = surf1.average_edge_length();
    double surf2_ael = surf2.average_edge_length();
   
    std::map< vertex_descriptor, double> map1,map2;
    
    vertex_vector surface1_vertices  = surf1.get_vertices_inside(surf2);  
    vertex_vector surface2_vertices  = surf2.get_vertices_inside(surf1);

    surface1_vertices  = surf1.get_vertices_outside(other,surface1_vertices);
    surface2_vertices  = surf2.get_vertices_outside(other,surface2_vertices);
    

    int iter =0;
    
    while (  !surface1_vertices.empty() or  !surface2_vertices.empty() )
    {
        surf1.get_adjacent_vertices(surface1_vertices);
        surf2.get_adjacent_vertices(surface2_vertices);
       
        map1 = surf1.get_shortest_edge_map(surface1_vertices,edge_movement);
        map2 = surf2.get_shortest_edge_map(surface2_vertices,edge_movement);

        surf1.adjust_vertices_in_region(map1.begin() ,map1.end());
        surf2.adjust_vertices_in_region(map2.begin() ,map2.end());
          
        surf1.smooth_laplacian_region(surface1_vertices.begin(), surface1_vertices.end(),smoothing);   
        surf2.smooth_laplacian_region(surface2_vertices.begin(), surface2_vertices.end(),smoothing);   
    
        //surface1_vertices = surf1.vertices_inside(surf2,surface1_vertices);
        //surface2_vertices = surf2.vertices_inside(surf1,surface2_vertices);
       //TODO
        if (surf1.num_self_intersections()>0)
        {
            surf1.collapse_edges(surf1_ael);
            surface1_vertices  = surf1.get_vertices_inside(surf2);
        
        }
        else 
            surface1_vertices  = surf1.get_vertices_inside(surf2,surface1_vertices);
        
        if (surf2.num_self_intersections()>0)
        {
            surf2.collapse_edges(surf2_ael);
            surface2_vertices  = surf2.get_vertices_inside(surf1);
        }
        else 
            surface2_vertices  = surf2.get_vertices_inside(surf1,surface2_vertices);
    


    
        surface1_vertices  = surf1.get_vertices_outside(other,surface1_vertices);
        surface2_vertices  = surf2.get_vertices_outside(other,surface2_vertices);

        if ( surf1.num_self_intersections()+surf2.num_self_intersections() >0 )
        {                
           std::cout << "Detected "<< surf1.num_self_intersections()+surf2.num_self_intersections() <<" self-intersections" << std::endl;
           std::cout << "Recommend use istropic_remeshing before continuation"  << std::endl;
           return false;
        }
        if ( iter++>max_iter)
        {
           std::cout << "Failed to converge in "<< max_iter <<" steps, terminating" << std::endl;
           return false;       
        }
    }
   return true ;
}

/**
 * @brief Separates two close surfaces outside a thrid surface.
 *
 * Separates the points of two surfaces that are outside a third surface iteratively by moving the vertices in the negative normal 
 * direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each 
 * vertex.This continues until the closest vertex is adjacent.
 *
 * @param surf1 first SVMTK surface class object  
 * @param surf2 second SVMTK surface class object 
 * @param other third SVMTK surface class obejct, the algorithm is not applied to vertices inside this surface.
 * @param double edge_movement the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param double smoothing the input of the laplacian smoothing that is applied after each iteration of vertex movement 
 * @param max_iter maximum number of iteration.
 * @return Surface class union of the modified 
 */
template< typename Surface> 
bool separate_close_surfaces(Surface& surf1 , Surface& surf2 , Surface& other, double edge_movement=-0.25, double smoothing=0.3, int max_iter=400)
{
   typedef typename Surface::vertex_vector vertex_vector;
   typedef typename Surface::vertex_descriptor vertex_descriptor;

   double surf1_ael = surf1.average_edge_length();
   double surf2_ael = surf2.average_edge_length();

   vertex_vector close_vertices_surf1,close_vertices_surf2;
   std::map< vertex_descriptor, double> map1,map2;

   // FIXME 
   close_vertices_surf1 = surf1.get_close_vertices(surf2); //get close vertices ??? 
   close_vertices_surf2 = surf2.get_close_vertices(surf1);

   
   int iter=0;
   while (  !close_vertices_surf1.empty() or  !close_vertices_surf2.empty() )
   {
   
        surf1.get_adjacent_vertices(close_vertices_surf1);
        surf2.get_adjacent_vertices(close_vertices_surf2);   
        
        close_vertices_surf1 = surf1.get_vertices_outside(other, close_vertices_surf1);
        close_vertices_surf2 = surf2.get_vertices_outside(other, close_vertices_surf2);

        map1 = surf1.get_shortest_edge_map(close_vertices_surf1,edge_movement);
        map2 = surf2.get_shortest_edge_map(close_vertices_surf2,edge_movement);

        surf1.adjust_vertices_in_region(map1.begin() ,map1.end());
        surf2.adjust_vertices_in_region(map2.begin() ,map2.end());

        surf1.smooth_laplacian_region( close_vertices_surf1.begin(),close_vertices_surf1.end(),smoothing);   
        surf2.smooth_laplacian_region( close_vertices_surf2.begin(),close_vertices_surf2.end(),smoothing);   
 
        //TODO
        if (surf1.num_self_intersections()>0)
        {
            surf1.collapse_edges(surf1_ael);
            close_vertices_surf1 = surf1.get_close_vertices(surf2);
        
        }
        else 
           close_vertices_surf1 = surf1.get_close_vertices(surf2,close_vertices_surf1);
        
        if (surf2.num_self_intersections()>0)
        {
            surf2.collapse_edges(surf2_ael);
            close_vertices_surf2 = surf2.get_close_vertices(surf1);     
        }
        else 
            close_vertices_surf2 = surf2.get_close_vertices(surf1,close_vertices_surf2);         
        



        //TODO
        //surf1.get_normal_vector_cluster( close_vertices_surf1);
        //surf2.get_normal_vector_cluster( close_vertices_surf2);
        
        if ( surf1.num_self_intersections()+surf2.num_self_intersections() >0 )
        {                
           std::cout << "Detected "<< surf1.num_self_intersections()+surf2.num_self_intersections() <<" self-intersections" << std::endl;
           std::cout << "Recommend use istropic_remeshing before continuation"  << std::endl;
           return false;
        }
        if ( iter++>max_iter)
        {
          
           std::cout << "Failed to converge in "<< max_iter <<" steps, terminating" << std::endl;
           return false;                
        }
   }
   return true;
}

/**
 * @brief Separates two close surfaces.
 * 
 * Separates two surfaces iteratively by moving the vertices in the negative normal 
 * direction that is determined by the multiplication of a negative value and the shortest edge length corresponding to each 
 * vertex. This continues until the closest vertex is adjacent. 
 *
 * @param surf1 first SVMTK surface class object  
 * @param surf2 second SVMTK surface class object 
 * @param double edge_movement the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param double smoothing the input of the laplacian smoothing that is applied after each iteration of vertex movement 
 * @param max_iter maximum number of iteration.
 * @return Surface class union of the modified 
 */
template< typename Surface> 
bool separate_close_surfaces(Surface& surf1 , Surface& surf2, double edge_movement=-0.25, double smoothing=0.1, int max_iter=400)
{
   typedef typename Surface::vertex_vector vertex_vector;
   typedef typename Surface::vertex_descriptor vertex_descriptor;

   double surf1_ael = surf1.average_edge_length(); // TODO: test
   double surf2_ael = surf2.average_edge_length();
   
   vertex_vector close_vertices_surf1,close_vertices_surf2;
   std::map< vertex_descriptor, double> map1,map2;

   close_vertices_surf1 = surf1.get_close_vertices(surf2);
   close_vertices_surf2 = surf2.get_close_vertices(surf1);


   // FIXME    
   surf1.get_adjacent_vertices(close_vertices_surf1);
   surf2.get_adjacent_vertices(close_vertices_surf2);
   int iter=0;
   while (  !close_vertices_surf1.empty() or  !close_vertices_surf2.empty() )
   {
        surf1.get_adjacent_vertices(close_vertices_surf1);
        surf2.get_adjacent_vertices(close_vertices_surf2);   
        
        map1 = surf1.get_shortest_edge_map(close_vertices_surf1,edge_movement);
        map2 = surf2.get_shortest_edge_map(close_vertices_surf2,edge_movement);

        surf1.adjust_vertices_in_region(map1.begin() ,map1.end());
        surf2.adjust_vertices_in_region(map2.begin() ,map2.end());

        surf1.smooth_laplacian_region( close_vertices_surf1.begin(),close_vertices_surf1.end(),smoothing);   
        surf2.smooth_laplacian_region( close_vertices_surf2.begin(),close_vertices_surf2.end(),smoothing);  
        

        if (surf1.num_self_intersections()>0)
        {
            surf1.collapse_edges(surf1_ael);
            close_vertices_surf1 = surf1.get_close_vertices(surf2); 
        }
        else 
           close_vertices_surf1 = surf1.get_close_vertices(surf2,close_vertices_surf1);
        
        if (surf2.num_self_intersections()>0)
        {
            surf2.collapse_edges(surf2_ael);
            close_vertices_surf2 = surf2.get_close_vertices(surf1);     
        }
        else 
            close_vertices_surf2 = surf2.get_close_vertices(surf1,close_vertices_surf2); 
            
        surf1.get_adjacent_vertices(close_vertices_surf1);
        surf2.get_adjacent_vertices(close_vertices_surf2);    
        //TODO WHY
        //surf1.get_normal_vector_cluster(close_vertices_surf1);
        //surf2.get_normal_vector_cluster(close_vertices_surf2);
        
        if ( surf1.num_self_intersections()+surf2.num_self_intersections() >0 )
        {                
           std::cout << "Detected "<< surf1.num_self_intersections()+surf2.num_self_intersections() <<" self-intersections" << std::endl;
           std::cout << "Recommend use istropic_remeshing before continuation"  << std::endl;
           return false;
        }

        if ( iter++>max_iter)
        {
           std::cout << "Failed to converge in "<< max_iter <<" steps, terminating" << std::endl;
           return false;                
        }
   }
   return true;
}

/**
 * @brief Takes surface union of surfaces that partially overlapp each other. 
 * 
 * Requires that the surfaces does overlap. 
 * Surface vertices inside the other surface are obtained. Then adjacent vertices
 * within a threshold angle of the vertex normal are added to the vertex vector.
 * The points corresponding to the vertex vector are iteratively moved in the 
 * vertex normal direction until there is suffient overlapp between surfaces.
 *        
 *
 * @tparam SVMTK Surface class.
 * @param Surface class 
 * @param Surface class
 * @param double clusterth lower bound of the cos angle between normal vectors that  
 * @param double edge_movement the multipler that with the smallest edge of a vertex that indicates the longest movement of that vertex.
 * @param int smoothing the number of taubin iterations after each iteration of vertex movement 
 * @param max_iter maximum number of iteration. 
 * @return Surface class union of the modified 
 */
template< typename Surface> 
std::shared_ptr<Surface> union_partially_overlapping_surfaces( Surface& surf1 , Surface& surf2, double clusterth, double edge_movement, int smoothing, int max_iter )
{

     // TODO add Precondition
     typedef typename Surface::vertex_vector vertex_vector;
     typedef typename Surface::vertex_descriptor vertex_descriptor;

     std::map< vertex_descriptor, double> map1,map2;
   
     vertex_vector surface1_vertices,surface2_vertices;

     surface1_vertices  = surf1.get_vertices_inside(surf2); 
     surface2_vertices  = surf2.get_vertices_inside(surf1);

     surf1.get_normal_vector_cluster(surface1_vertices,clusterth);
     surf2.get_normal_vector_cluster(surface2_vertices,clusterth);
     surf1.get_adjacent_vertices(surface1_vertices);
     surf2.get_adjacent_vertices(surface2_vertices);
     int iter =0;
    
     while (  !surface1_vertices.empty() or  !surface2_vertices.empty() )   // 
     {          
           map1 = surf1.get_shortest_edge_map(surface1_vertices,edge_movement);
           map2 = surf2.get_shortest_edge_map(surface2_vertices,edge_movement);

           surf1.adjust_vertices_in_region(map1.begin() ,map1.end());
           surf2.adjust_vertices_in_region(map2.begin() ,map2.end());
           
           surf1.smooth_taubin_region(surface1_vertices.begin(), surface1_vertices.end(),smoothing);   
           surf2.smooth_taubin_region(surface2_vertices.begin(), surface2_vertices.end(),smoothing); 

           surface1_vertices  = surf1.get_vertices_outside(surf2,surface1_vertices); 
           surface2_vertices  = surf2.get_vertices_outside(surf1,surface2_vertices);
           if ( iter++>max_iter)
           {
              break;
           }
     }
     std::shared_ptr<Surface> result(new Surface()); 
     CGAL::Polygon_mesh_processing::corefine_and_compute_union(surf1.get_mesh(), surf2.get_mesh(), result->get_mesh());     
     return result;
}

/**
 * \class Surface 
 *
 * SVMTK Surface class is used to handle operations related to triangualted 
 * surfaces in 3D CGAL Surface Mesh 
 * @see ()[]
 *   
 * 
 * It should be noted that most operations can also be used on CGAL
 * Polyhedral Surface 
 * @see()[]
 * 
 *                       -     
 *
 *
 */ 
class Surface
{
   public:
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel; 

    typedef Kernel::Point_3 Point_3;
    typedef Kernel::Plane_3 Plane_3;
    typedef Kernel::FT FT; 
    typedef Kernel::Sphere_3 Sphere_3; 
    typedef Kernel::Vector_3 Vector_3;
    typedef Kernel::Triangle_3 Triangle_3;

    typedef CGAL::Surface_mesh<Point_3> Mesh;
    typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
    typedef CGAL::Side_of_triangle_mesh<Mesh,Kernel> Inside; 
    typedef CGAL::Polyhedron_3<Kernel> Polyhedron;   
    typedef CGAL::Vertex_around_target_circulator<Mesh> HV_const_circulator; 


    typedef boost::graph_traits<Mesh>::edge_descriptor             edge_descriptor;
    typedef boost::graph_traits<Mesh>::halfedge_descriptor         halfedge_descriptor;
    typedef boost::graph_traits<Mesh>::face_descriptor             face_descriptor;
    typedef boost::graph_traits<Mesh>::vertex_descriptor           vertex_descriptor; 
    typedef boost::graph_traits<Mesh>::vertices_size_type          size_type;
    typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type   Vertex_point_pmap; 

    typedef std::vector<Point_3>               Polyline;
    typedef std::vector<Polyline>              Polylines;
    typedef std::vector<std::size_t>           Face;
    typedef std::vector<Surface::Point_3>      point_vector;
    typedef std::vector<vertex_descriptor>     vertex_vector;       
    typedef std::vector<face_descriptor>       face_vector; 

    typedef Mesh::Vertex_index                 Index;
    
    //TODO: Rename?
    typedef std::map<vertex_descriptor,Vector_3> vertex_vector_map;
    typedef std::map<vertex_descriptor,double> vertex_scalar_map;     
       
    typedef CGAL::Search_traits_3<Kernel>                                                Traits_base;
    typedef CGAL::Search_traits_adapter<vertex_descriptor,Vertex_point_pmap,Traits_base> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>                                   K_neighbor_search;
    typedef K_neighbor_search::Tree                                                      Tree;
    typedef Tree::Splitter                                                               Splitter;
    typedef K_neighbor_search::Distance                                                  Distance;

    /* -- Constructors -- */

    Surface(){} 
    Surface(Polyhedron &polyhedron); 
    Surface(std::vector<Point_3>& points,std::vector<Face>& faces ); 
    Surface(const std::string  filename);
    Surface(const Surface &other) {this->mesh=other.mesh; } 
    Surface(const std::shared_ptr<Surface> surf) {this->mesh=surf.get()->get_mesh(); } 
    ~Surface(){}
    
    Surface &operator=(Surface &other){ this->mesh= other.mesh; return *this; }
    bool is_point_inside(Point_3 point_3);


    void make_cone_(    double x0, double y0, double  z0,  double x1, double y1, double z1, double r0, double edge_length) ;
    void make_cylinder( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0, double edge_length) ;
    void make_cone(     double x0, double y0, double  z0,  double x1, double y1, double z1, double r0,double r1, double edge_length) ; 
    void make_cube(     double x0, double y0, double  z0,  double x1, double y1, double z1, double edge_length);
    void make_sphere(   double x0, double y0, double  z0,  double r0, double edge_length); 
    void make_circle_in_plane(Point_3 point, Vector_3 vector, double radius, double edge_length);
    
   
    void make_circle_in_plane(double px, double py, double pz, double vx, double vy, double vz, double radius, double edge_length)
    {
         make_circle_in_plane(Point_3(px,py,pz),Vector_3(vx,vy,vz), radius, edge_length);
    }

    void make_cube(Point_3 p0, Point_3 p1, double edge_length )
    {
         make_cube(CGAL::to_double(p0.x()), CGAL::to_double(p0.y()),CGAL::to_double(p0.z()),
                   CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),CGAL::to_double(p1.z()),
                   edge_length);
    }
     
    void make_cylinder(Point_3 p0,Point_3 p1, double r0, double edge_length)
    {
         make_cylinder(CGAL::to_double(p0.x()), CGAL::to_double(p0.y()),CGAL::to_double(p0.z()),
                       CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),CGAL::to_double(p1.z()),
                       r0,edge_length);
    }
    void make_cone(Point_3 p0,Point_3 p1, double r0, double r1, double edge_length)
    {
         make_cone(CGAL::to_double(p0.x()), CGAL::to_double(p0.y()),CGAL::to_double(p0.z()),
                   CGAL::to_double(p1.x()), CGAL::to_double(p1.y()),CGAL::to_double(p1.z()),
                   r0,r1,edge_length);
    }
    void make_sphere(Point_3 p0, double r0, double edge_length)
    {
         make_sphere(CGAL::to_double(p0.x()), CGAL::to_double(p0.y()),CGAL::to_double(p0.z()), 
                     r0, edge_length);
    }
   
    bool surface_intersection(Surface other);
    bool surface_difference(Surface other);
    bool surface_union(Surface other);

    Mesh& get_mesh() {return mesh;}
    void clear(){ mesh.clear();}

    int num_faces()    const {return mesh.number_of_faces();}
    int num_edges()    const {return mesh.number_of_edges();}
    int num_vertices() const {return mesh.number_of_vertices();}
    int num_self_intersections();

    void save(const std::string outpath);
   
    bool is_empty(){ return mesh.is_empty();}
   
    template< typename Polyhedron_3>  
    void get_polyhedron(Polyhedron_3 &polyhedron_3 ){ assert_non_empty_mesh(); CGAL::copy_face_graph(mesh,polyhedron_3);}    

    void get_normal_vector_cluster( vertex_vector &vertices,double angle_in_degree=36.87);
    void smooth_shape(double time,int nb_iterations);

    int  fill_holes();                       
    bool triangulate_faces();                 

    void split_edges(double  target_edge_length);
    int  collapse_edges(const double target_edge_length);
    int  collapse_edges();    

    void isotropic_remeshing(double target_edge_length, unsigned int nb_iter, bool protect_border);

    void adjust_boundary(const double c);
    void smooth_laplacian(const double c, int iter);
    void smooth_taubin(const size_t nb_iter); 


    bool clip(double a,double b, double c ,double d, bool preserve_manifold);
    bool clip(Point_3 point, Vector_3 vector , bool preserve_manifold);
    bool clip(Plane_3 plane ,bool preserve_manifold);
    bool clip(Surface other,bool invert,bool preserve_manifold);
    bool clip(Point_3 point, Vector_3 vector, double radius, bool invert  ,bool preserve_manifold);
        
    std::vector<std::pair<Point_3, Vector_3>> get_points_with_normal();

    // TODO rename mesh_slice -> get_slice
    template< typename Slice>
    std::shared_ptr<Slice> mesh_slice(double x1,double x2, double x3 ,double x4) ;
    template<typename Slice>
    std::shared_ptr<Slice> mesh_slice(Plane_3 plane) ;

    std::pair<double,double> span(int direction);

    std::shared_ptr< Surface > convex_hull();

    // tempalte used for Vetex::range also 
    template<typename InputIterator>   
    void adjust_vertices_in_region(InputIterator begin, InputIterator end, const double c);
    template<typename InputIterator>
    void smooth_laplacian_region(InputIterator begin, InputIterator end  ,const double c);
    template<typename InputIterator>
    void smooth_taubin_region(InputIterator begin, InputIterator end     ,const size_t iter);

    vertex_scalar_map get_shortest_edge_map(const Surface::vertex_vector &vector,const double adjustment =-0.5) ;

    void adjust_vertices_in_region(vertex_scalar_map::iterator begin, vertex_scalar_map::iterator end );
    void adjust_vertices_in_region(vertex_vector_map::iterator begin, vertex_vector_map::iterator end );

    void smooth_laplacian_region(vertex_vector_map::iterator begin, vertex_vector_map::iterator end,const double c);
    
    vertex_vector get_vertices();
    point_vector  get_points(); 
    point_vector  get_points(Surface::vertex_vector &vertices); 

    template<int A=0>
    vertex_vector_map get_close_vertices_with_direction(Surface& other, double adjustment);

    Polyline get_shortest_surface_path(Point_3 source, Point_3 target);
    Polyline get_shortest_surface_path(double x0, double y0, double z0, double x1, double y1, double z1);

    void reconstruct(double angular_bound=20, double radius_bound=0.1, double distance_bound=0.1);

    Polylines mean_curvature_flow();
    
    
    Polylines polylines_in_plane(Plane_3 plane_3);
    
    Polylines polylines_in_plane(Point_3 p, Vector_3 v)
    {
        return polylines_in_plane(Plane_3(p,v));
    }
    
    
    
    template<typename Implicit_function>  
    void implicit_surface(Implicit_function implicit_function,
             double bounding_sphere_radius,
             double angular_bound=30.,
             double radius_bound=0.1 ,
             double distance_bound=0.1   );

    void set_outward_face_orientation();
 
    std::shared_ptr<Surface> cylindric_extension(const Point_3& p1,double radius, double length, double edge_length ,bool normal=true );      
    std::shared_ptr<Surface> cylindric_connection(Surface other, double radius, double edge_length);              
                   
    vertex_vector get_closest_vertices(Point_3 p1, int num = 8);
    point_vector  get_closest_points(Point_3 p1, int num=8);
 
    int keep_largest_connected_component();
    
    Point_3 centeroid(){return CGAL::Polygon_mesh_processing::centroid(mesh);}   
    double volume(){return CGAL::to_double(CGAL::Polygon_mesh_processing::volume(mesh));}
    double area(){return CGAL::to_double(CGAL::Polygon_mesh_processing::area(mesh));}
    double distance_to_point(Surface::Point_3 point);
    std::string CGAL_precondition_evaluation(Surface other);
    bool does_bound_a_volume();
   
    std::vector<std::pair<Triangle_3, std::pair<int,int>>>  surface_segmentation(int nb_of_patch_plus_one=1,double angle_in_degree=85);
 
    // TODO : templates or "overloaded functions" same funtion one line difference ? 
    template <int A=0> 
    vertex_vector get_close_vertices(Surface &other);
    
    template <int A=0> 
    vertex_vector get_close_vertices(Surface &other,vertex_vector &mvertices);
    vertex_vector get_narrow_gaps();
    
    template<CGAL::Bounded_side A , CGAL::Bounded_side B> 
    std::pair<bool,int>  manipulate_vertex_selection(Surface &other, double adjustment,  double smoothing, int max_iter);
  
    template< int A=0>
    std::pair<bool,int>  manipulate_vertex_selection(Surface &other, double adjustment,  double smoothing, int max_iter);  
    
    template<int A=0>
    std::pair<bool,int> manipulate_vertex_selection_with_direction(Surface &other, double adjustment, int max_iter);
    
    
    
    template< CGAL::Bounded_side A , CGAL::Bounded_side B>
    vertex_vector get_vertices_with_property(Surface &other) ; 
        
    template< CGAL::Bounded_side A , CGAL::Bounded_side B>
    vertex_vector get_vertices_with_property(Surface &other, vertex_vector &vertices) ;

    /* --  --*/
    vertex_vector get_vertices_outside(Surface &other); 
    vertex_vector get_vertices_inside(Surface &other); 
    vertex_vector get_vertices_outside(Surface &other    , vertex_vector &vertices ); 
    vertex_vector get_vertices_inside(Surface &other     , vertex_vector &vertices  ); 
    vertex_vector get_vertices_on_boundary(Surface &other, vertex_vector &vertices );
    
    // TODO overload instead 
    void get_adjacent_vertices(vertex_vector &vertices); 
    void get_adjacent_vertices_with_direction(vertex_vector_map &vertices_wdir);
    
    std::pair<bool,int> embed(   Surface& other, double adjustment=-0.5, double smoothing=0.3, int max_iter=400); 
    
    std::pair<bool,int> enclose( Surface& other, double adjustment=0.5 , double smoothing=-0.3, int max_iter=400); 
    
    // TODO: REMOVE
    std::pair<bool,int> separate(Surface& other, double adjustment=0.5,int max_iter=400); 

    std::pair<bool,int>  separate_narrow_gaps(double adjustment, double smoothing, int max_iter) ;
    std::pair<bool,int>  separate_close_vertices(double adjustment, int max_iter) ;

    /**
     * @brief Throws error if the class variable mesh is empty
     * 
     * @param none 
     * @return true if mesh is non empty, otherwise false 
     */
    bool assert_non_empty_mesh()
    { 
      if (is_empty())
        throw EmptyMeshError("Surface is empty") ;
      else 
        return true;
      return false;
    }  
   
   
    template< CGAL::Bounded_side A , CGAL::Bounded_side B>
    std::pair<bool,bool> check_vertices(Surface &other);
     
    //void  fix_self_intersections();


    /*void triangulate_hole()
    {
         CGAL::Polygon_mesh_processing::triangulate_hole(mesh);	
    }*/
   /**
    * @brief Creates a cylinder that extends from the surface.   
    * 
    * @overload 
    */   
    std::shared_ptr<Surface> cylindric_extension(double x, double y ,double z, double radius, double length,double edge_length, bool normal)
                   { return cylindric_extension(Point_3(x,y,z),radius,length,edge_length,normal);}


   double average_edge_length();
   protected:
    Mesh mesh;
    

};

/**
 * TODO : implemnt CGAL function
 * Experimental 
 *
 * @param vertices a vector of vertices.
 * @return none 
 */
/*inline void Surface::fix_self_intersections()
{

}*/

/* -- Vertex Queries -- */

/**
 * Updates argument vertices to include connected adjacent vertices.
 *
 * @param vertices a vector of vertices.
 * @return none 
 */
inline void Surface::get_adjacent_vertices(Surface::vertex_vector &vertices) 
{
   
    vertex_vector temp;
    for ( vertex_descriptor vit : vertices )  
    {
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(vit), mesh), done(vbegin);
        do
        {
         if(std::find(vertices.begin(), vertices.end(), *vbegin) == vertices.end() ) 
             temp.push_back(*vbegin);
         *vbegin++;
        }while(vbegin!=done);
    }
    vertices.insert(vertices.end(), temp.begin(), temp.end());
}

/**
 * Updates argument vertices to include connected adjacent vertices.
 *
 * @param vertices a vector of vertices.
 * @return none 
 */
inline void Surface::get_adjacent_vertices_with_direction(Surface::vertex_vector_map &vertices_wdir) 
{

    vertex_vector_map temp;
    for ( auto vit : vertices_wdir )  
    {
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(vit.first), mesh), done(vbegin);
        do
        {
         if( vertices_wdir.find(*vbegin) == vertices_wdir.end()) 
            temp[*vbegin]=0.5*vit.second;
         *vbegin++;
        }while(vbegin!=done);
    }
    vertices_wdir.insert(temp.begin(), temp.end());
}

/**
 * @brief Returns the vertices of the surface mesh
 * @param none 
 * @return vector of surface vertices 
 */
inline Surface::vertex_vector Surface::get_vertices() 
{
   assert_non_empty_mesh();
   vertex_vector result;
   for ( vertex_descriptor vit : mesh.vertices()){result.push_back(vit);}
   return result;
   
}

/**
 * @brief Finds and return vertices that are inside another SVMTK Surface class.
 
 * @tparam CGAL::Bounded_side A
 * @tparam CGAL::Bounded_side B
 * @param other SVMTK Surface class 
 * @param vertices vector of this-> surface mesh vertices   
 * @return result vector of vertices 
 * @overload 
 */
template< CGAL::Bounded_side A , CGAL::Bounded_side B>
inline Surface::vertex_vector Surface::get_vertices_with_property(Surface &other)
{
   assert_non_empty_mesh();
   vertex_vector result;
   Surface::Inside is_inside_query(other.get_mesh()); 
   for ( vertex_descriptor vit : mesh.vertices() )
   {
      CGAL::Bounded_side res =  is_inside_query(mesh.point(vit));
      if (res == A or res == B)
         result.push_back(vit);
   }
   return result;
}

/**
 * @brief Finds and return vertices that are inside another SVMTK Surface class.

 * @tparam CGAL::Bounded_side A
 * @tparam CGAL::Bounded_side B
 * @param other SVMTK Surface class 
 * @param vertices vector of this-> surface mesh vertices   
 * @return result vector of vertices 
 * @overload 
 */
template< CGAL::Bounded_side A , CGAL::Bounded_side B>
inline Surface::vertex_vector Surface::get_vertices_with_property(Surface &other,Surface::vertex_vector &vertices)
{
   assert_non_empty_mesh();
   vertex_vector result;
   Surface::Inside is_inside_query(other.get_mesh()); 
   for ( vertex_descriptor vit : vertices )
   {
      CGAL::Bounded_side res =  is_inside_query(mesh.point(vit));
      if (res == A or res == B)
         result.push_back(vit);
   }
   return result;
}

/**
 * @brief checks vertices position property relative to surface.
 *
 * This function is used to investigate if a surface is partially 
 * inside and outside of another surface.
 *  
 * @tparam CGAL::Bounded_side CGAL Enum 
 * @see [CGAL::Bounded_side] (https://doc.cgal.org/latest/Kernel_23/group__kernel__enums.html)
 *
 * @param other SVMTK Surface class object 
 * @return pair of bool
 * 
 */ 
 template< CGAL::Bounded_side A , CGAL::Bounded_side B>
inline std::pair<bool,bool> Surface::check_vertices(Surface &other)
{
   assert_non_empty_mesh();
   bool  query1 = false; 
   bool  query2 = false; 
   Surface::Inside is_inside_query(other.get_mesh()); 
   for ( vertex_descriptor vit : mesh.vertices() )
   {
      CGAL::Bounded_side res =  is_inside_query(mesh.point(vit));
      if (res == A)                                
          query1 = true;
      if (res == B )
          query2 = true;
   }
   return std::make_pair(query1,query2);
}

/** 
 * @brief Finds and returns vertices that are outside another SVMTK Surface class.
 * @param other SVMTK Surface class
 * @param vertices vector of surface mesh vertices   
 * @return result vector of vertices 
 * @overload 
 */
inline Surface::vertex_vector Surface::get_vertices_outside(Surface &other,Surface::vertex_vector &vertices) 
{
   return this->get_vertices_with_property<CGAL::ON_UNBOUNDED_SIDE,
                              CGAL::ON_UNBOUNDED_SIDE>(other, vertices) ;
}

/** 
 * @brief Finds and returns vertices that are inside another SVMTK Surface class.
 * @param other SVMTK Surface class
 * @param vertices vector of surface mesh vertices   
 * @return result vector of vertices 
 * @overload 
 */
inline Surface::vertex_vector Surface::get_vertices_inside(Surface &other,Surface::vertex_vector &vertices) 
{
   return this->get_vertices_with_property<CGAL::ON_BOUNDED_SIDE,
                              CGAL::ON_BOUNDARY>(other, vertices);
}

/** 
 * @brief Finds and returns vertices that are on the boundary of another SVMTK Surface class.
 * @param other SVMTK Surface class
 * @param vertices vector of surface mesh vertices   
 * @return result vector of vertices 
 * @overload 
 */
inline Surface::vertex_vector Surface::get_vertices_on_boundary(Surface &other,Surface::vertex_vector &vertices) 
{
   return this->get_vertices_with_property<CGAL::ON_BOUNDARY,CGAL::ON_BOUNDARY>(other, vertices) ;
}

/** 
 * @brief Finds and returns vertices that are outside another SVMTK Surface class.
 * @param other SVMTK Surface class
 * @return result vector of vertices 
 * @overload 
 */
inline Surface::vertex_vector Surface::get_vertices_outside(Surface &other) 
{
   vertex_vector vertices = get_vertices();
   return this->get_vertices_with_property<CGAL::ON_UNBOUNDED_SIDE,
                              CGAL::ON_UNBOUNDED_SIDE>(other,vertices) ;
}

/** 
 * @brief Finds and returns vertices that are inside another SVMTK Surface class.
 * @param other SVMTK Surface class
 * @return result vector of vertices 
 * @overload 
 */
inline Surface::vertex_vector Surface::get_vertices_inside(Surface &other) 
{
   vertex_vector vertices = get_vertices();
   return this->get_vertices_with_property<CGAL::ON_BOUNDED_SIDE,
                              CGAL::ON_BOUNDARY>(other,vertices) ;
}

/**
 * @brief Finds and returns surface mesh vertices that are close to another SVMTK Surface Class
 *
 * @note removes all isolated vertices to avoid error in algorithm
 * @param other SVMTK Surface Class.
 * @return results a vector of vertices. 
 */
template <int A> 
inline Surface::vertex_vector Surface::get_close_vertices(Surface &other)
{
   assert_non_empty_mesh();
   
   CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);

   Surface::vertex_vector results;
   Vertex_point_pmap vppmap = get(CGAL::vertex_point,other.get_mesh());  

   Tree tree(vertices(other.get_mesh()).begin(),   // other.get_mesh().vertices().begin()?? 
            vertices(other.get_mesh()).end(),
            Splitter(),
            Traits(vppmap)
   );
   Distance tr_dist(vppmap);
   FT distance, edgeL;
   Point_3 closest;
   Vector_3 direction,normal;
   bool flag;

   for (vertex_descriptor vit : mesh.vertices())
   {
        flag = true;
        
        K_neighbor_search search(tree, mesh.point(vit), 2,0,true,tr_dist); 
        
        closest = other.get_mesh().point((search.begin()+A)->first); 

        Point_3 current = mesh.point(vit);
        
        direction =  Vector_3(closest,current);
 
        normal =  CGAL::Polygon_mesh_processing::compute_vertex_normal(vit,mesh);
     
        distance = direction.squared_length();
        //direction= direction/distance; 
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(vit), mesh), done(vbegin);

        if (normal*direction<=0) // FIXME
        {
           continue;
        }

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

        if (flag){ results.push_back(vit);} // both vertices are affected
   }
   return results;
}

/**
 * @brief Finds and returns surface mesh vertices that are close to another SVMTK Surface Class
 *
 * @note removes all isolated vertices to avoid error in algorithm
 * @param other SVMTK Surface Class.
 * @return results a vector of vertices. 
 */
template <int A> 
inline Surface::vertex_vector Surface::get_close_vertices(Surface &other, Surface::vertex_vector &mvertices)
{
   assert_non_empty_mesh();
   
   CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);

   Surface::vertex_vector results;
   Vertex_point_pmap vppmap = get(CGAL::vertex_point,other.get_mesh());  

   Tree tree(vertices(other.get_mesh()).begin(),
            vertices(other.get_mesh()).end(),
            Splitter(),
            Traits(vppmap)
   );
   Distance tr_dist(vppmap);
   FT distance, edgeL;
   Point_3 closest;
   Vector_3 direction,normal;
   bool flag;

   for (vertex_descriptor vit : mvertices )
   {
        flag = true;
        
        K_neighbor_search search(tree, mesh.point(vit), 2,0,true,tr_dist); 
        
        closest = other.get_mesh().point((search.begin()+A)->first); 

        Point_3 current = mesh.point(vit);
        
        direction =  Vector_3(closest,current);
 
        normal =  CGAL::Polygon_mesh_processing::compute_vertex_normal(vit,mesh);
     
        distance = direction.squared_length();

        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(vit), mesh), done(vbegin);

        if (normal*direction<=0) //FIXME diables self intersections 
        {
           continue;
        }

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

        if (flag){ results.push_back(vit);} // both vertices are affected
   }
   return results;
}

/**
 * @brief Finds close non-adjacent vertices and returns a map with an direction to move them apart.
 *
 * Returns vertices that are closer to another SVMTK surface object than any of the edge connected adjacent vertices.
 * @note removes all isolated vertices to avoid error in algorithm
 * @param other a SVMTK Surface class object.
 * @param adjustment a negative multiplier for the direction between close vertices
 * @return results a std::map with vertices as keys and Vector_3 as value.  
 *  
 */
 template < int A>
inline Surface::vertex_vector_map Surface::get_close_vertices_with_direction(Surface& other, double adjustment)
{  
   //if (adjustment>=0) 
   //   throw InvalidArgumentError("The adjust must be negative, i.e. contraction.");

   assert_non_empty_mesh();
   CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh); 
     
   vertex_vector_map results;
   
   Vertex_point_pmap vppmap = get(CGAL::vertex_point,other.get_mesh());
   Tree tree(vertices(other.get_mesh()).begin(),
            vertices(other.get_mesh()).end(),
            Splitter(),
            Traits(vppmap)
   );

   Distance tr_dist(vppmap);

   FT distance,edgeL;
   Point_3 closest;
   Vector_3 direction,normal;
   bool flag;
   for (vertex_descriptor vit : mesh.vertices())
   {
        flag = true;
    
        K_neighbor_search search(tree, mesh.point(vit), 2,0,true,tr_dist); 
        
        closest = other.get_mesh().point((search.begin()+A)->first); 

        Point_3 current = mesh.point(vit);

        direction =  Vector_3(closest,current);

     
        distance = direction.squared_length();

        if (distance==0) 
           continue;
           
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(vit),mesh), done(vbegin);
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
    
        if (flag)
        {
           results[vit] = adjustment*direction; 
        }
   }
   return  results;   
}

/** 
 * @brief Separates non-adjacent surface mesh vertices.
 * 
 * Separates non-adjacent surface mesh vertices so that the nearest surface mesh vertices are adjacent. 
 * @note1 In case of inside-out surfaces with self-intersections, the algorithm may fail to produce  For the algorithm to have optimal effect.
 * @note2 removes all isolated vertices to avoid error in algorithm
 * 
 * @param adjustment multiplier of the edge movement.
 * @return number of adjusted vertices.  
 */
inline Surface::vertex_vector Surface::get_narrow_gaps() 
{
   assert_non_empty_mesh();     
   return get_close_vertices<1>(*this);
}

/* --- Vertex Algorithm -- */

/**
 * @brief Moves the vertex position of a selection of Surface vertices dependent
 * on different queries, like inside or outside another surface.
 *
 * @tparam CGAL::Bounded_side
 * @see Surface::get_vertices_ 
 *
 * @param other SVMTK Surface class object
 * @param adjustment times the shortest connected edge length gives the limit to the vertex movement. 
 * @param smoothing laplacian smoothing parameter for smoothing after each iteration. 
 * @param max_iter the maximum number of iteration for the while. 
 *
 * @return std::pair<bool,int>. True if algorithm complets before max_iter is 
 * reached. The integer is the number of vertex manipulation that
 * is completed.  
 *
 * FIXME
 */
template< CGAL::Bounded_side A , CGAL::Bounded_side B>
inline std::pair<bool,int> Surface::manipulate_vertex_selection(Surface &other, double adjustment,  double smoothing, int max_iter)
{
   std::map< vertex_descriptor, double> se_map;
   vertex_vector vertices = get_vertices_with_property<A,B>(other) ;
   int after=0;
   int before = vertices.size();
      
   int iter = 0;
   while ( !vertices.empty())
   {     
   

       se_map = get_shortest_edge_map(vertices,adjustment);  
         
       adjust_vertices_in_region(se_map.begin() ,se_map.end());
       
       smooth_laplacian_region(vertices.begin(),vertices.end(),smoothing);  
            
       vertices = this->get_vertices_with_property<A,B>(other);

       after = vertices.size();    

       if (++iter>max_iter and after>0)
          return std::make_pair(false, after);
       if ( after > before ) 
          throw AlgorithmError("The aglorithm failed, select a smaller adjustment parameter"); 
   }
   return std::make_pair(true, after);
}

/**
 * FIXME 
 * @brief Moves the vertex position of a selection of Surface vertices dependent
 * on adjaceny to other vertices, itself and other surfaces,
 *
 * @tparam int indicating the starting point of the adjacency search. If other 
 * is (*this) then A=1 so that vertex does not choose itself as a 
 * close vertex.  
 * @see Surface::get_close_vertices
 * 
 * @note The algorithm will collapse edges that are below the threshold of 
 * 0.4 of the average edge size.
 *
 * @param other SVMTK Surface class object
 * @param adjustment times the shortest connected edge length gives the limit to the vertex movement. 
 * @param smoothing laplacian smoothing parameter for smoothing after each iteration. 
 * @param max_iter the maximum number of iteration for the while. 
 *
 * @return std::pair<bool,int>. True if algorithm complets before max_iter is 
 * reached. The integer is the number of vertex manipulation that
 * is completed.  
 * 
 */
template<int A>
inline std::pair<bool,int> Surface::manipulate_vertex_selection(Surface &other, double adjustment,  double smoothing, int max_iter)
{


   std::map< vertex_descriptor, double> se_map;

   auto vertices = get_close_vertices<A>(other);
   
   
   //get_adjacent_vertices(vertices);  
        
   int after=0;
   int before = vertices.size();
   int iter = 0;
   while ( !vertices.empty())
   {         

       se_map = get_shortest_edge_map(vertices,adjustment);  
         
       adjust_vertices_in_region(se_map.begin() ,se_map.end());
  
       smooth_laplacian_region(vertices.begin(),vertices.end(),smoothing);      //TODO ReMOVE  ADDeD 
       
       vertices = get_close_vertices<A>(other,vertices);

       //get_adjacent_vertices(vertices);        
        
       after = vertices.size();  
            
       if (++iter>max_iter and after>0)
          return std::make_pair(false, after);
       if ( after > before ) 
          throw AlgorithmError("The aglorithm failed, select a smaller adjustment parameter"); 
   }
   return std::make_pair(true, after);
}

/**
 * @brief Moves the vertex position of a selection of Surface vertices dependent
 * on adjaceny to other vertices, itself and other surfaces,
 *
 * @tparam int indicating the starting point of the adjacency search. If other 
 * is (*this) then A=1 so that vertex, itself, is not chosen as the 
 * closest vertex.  
 * @see Surface::get_close_vertices
 * 
 * @param other SVMTK Surface class object
 * @param adjustment times the shortest connected edge length gives the limit to the vertex movement. 
 * @param smoothing laplacian smoothing parameter for smoothing after each iteration. 
 * @param max_iter the maximum number of iteration for the while. 
 *
 * @return std::pair<bool,int>. True if algorithm complets before max_iter is 
 * reached. The integer is the number of vertex manipulation that
 * is completed.  
 * 
 */
template<int A>
inline std::pair<bool,int> Surface::manipulate_vertex_selection_with_direction(Surface &other, double adjustment, int max_iter)
{


   vertex_vector_map vertices = get_close_vertices_with_direction<A>(other,adjustment);
   
   int after=0;
   int before = vertices.size();
   
   int iter = 0;
   while ( !vertices.empty())
   {         
       get_adjacent_vertices_with_direction(vertices);       //FIXME

       adjust_vertices_in_region(vertices.begin() ,vertices.end());

       smooth_laplacian_region(vertices.begin(),vertices.end(),0.5*adjustment);      //TODO ReMOVE  

       vertices = get_close_vertices_with_direction<A>(other,adjustment); //FIXME ERROR

       after = vertices.size();  
              
       if (++iter>max_iter and after>0)
          return std::make_pair(false, after);
       if ( after > before ) 
          throw AlgorithmError("The aglorithm failed, select a smaller adjustment parameter"); 
   }
   return std::make_pair(true, after);
}

/** 
 * @brief Separates non-adjacent surface mesh vertices in the negative surface normal direction, i.e.
 * contraction.
 * 
 * Separates non-adjacent surface mesh vertices in the negative direction of the vertex normal,
 * so that the nearest surface mesh vertices are adjacent. 
 * Self-intersection can cause the algorithm to fail. 
 * All isolated vertices are removed to avoid error.
 *
 * @note There are several required factors to finish before the maximum iteration is exceeded. However,
 * only a few iteration may be nesscary to obtain good result. 
 * 
 * @param adjustment multiplier of the edge movement.
 * @return std::pair<bool,int> completion before maximum iteration reached, and
 * number of adjusted vertices.  
 */
inline std::pair<bool,int>  Surface::separate_narrow_gaps(double adjustment, double smoothing, int max_iter) 
{
   assert_non_empty_mesh();
   CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh); 
   
   if (adjustment >0)
      throw  InvalidArgumentError("Adjusment must be negative.");
   return manipulate_vertex_selection<1>(*this,adjustment,smoothing,max_iter);
}

/** 
 * @brief Separates non-adjacent surface mesh vertices.
 * 
 * Separates non-adjacent surface mesh vertices so that the nearest surface mesh vertices are adjacent. 
 * 
 * Self-intersection can cause the algorithm to fail. 
 * All isolated vertices are removed to avoid error.
 *
 * @note There are several required factors to finish before the maximum iteration is exceeded. However,
 * only a few iteration may be nesscary to obtain good result. 
 * 
 * @param adjustment multiplier of the edge movement.
 * @return std::pair<bool,int> completion before maximum iteration reached, and
 * number of adjusted vertices.  
 */
inline std::pair<bool,int>  Surface::separate_close_vertices(double adjustment, int max_iter) 
{
   assert_non_empty_mesh();
   CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh); 
   
   if (adjustment <0)
      throw  InvalidArgumentError("Adjusment must be positive.");
   return manipulate_vertex_selection_with_direction<1>(*this,adjustment,max_iter);
}

/**
 * @brief Moves all vertices so that they are inside another surface. 
 *
 * @param other SVMTK Surface class object. FIXME
 * @param adjustment 
 * @param smoothing
 * @param max_iter 
 *
 * @return std::pair<bool,int> completion before maximum iteration reached, and
 * number of adjusted vertices. 
 *
 * @throws InvalidArgumentError if adjustment is positive,i.e.
 * expansion instead of contraction.
 */

inline std::pair<bool,int> Surface::embed(Surface& other, double adjustment,  double smoothing, int max_iter) 
{
   if (adjustment >=0)
      throw  InvalidArgumentError("Adjusment must be negative.");
      
   return manipulate_vertex_selection<CGAL::ON_UNBOUNDED_SIDE,CGAL::ON_BOUNDARY>(other, adjustment,smoothing,max_iter);
   
}

/**
 * @brief Moves all vertices so that they are outside another surface.  
 * Precondition that the surfaces intersect. FIXME
 *
 * @param other SVMTK Surface class object.
 * @param adjustment 
 * @param smoothing
 * @param max_iter
 *
 * @return std::pair<bool,int> completion before maximum iteration reached, and
 * number of adjusted vertices. 
 *
 * @throws InvalidArgumentError if adjustment is negative,i.e.
 *  contraction instead of expansion.
 * 
 */
inline std::pair<bool,int> Surface::enclose(Surface& other, double adjustment,double smoothing, int max_iter) 
{
     
   if (adjustment <=0)
      throw  InvalidArgumentError("Adjusment must be poisitive.");
   return manipulate_vertex_selection<CGAL::ON_BOUNDED_SIDE,CGAL::ON_BOUNDARY>(other, adjustment,smoothing,max_iter);
   
}

/**
 * @brief  Moves all vertices that are close to another surface in opposite direction.
 * Precondition that the other surface is either enclosed or embedded in the class
 * surface. 
 *
 * @param other SVMTK Surface class object.
 * @param adjustment 
 * @param smoothing
 * @param max_iter
 *
 
 * @return std::pair<bool,int> completion before maximum iteration reached, and
 * number of adjusted vertices. 
 */
inline std::pair<bool,int> Surface::separate(Surface& other, double adjustment, int max_iter)
{ 

    std::pair<bool,bool> queries = check_vertices<CGAL::ON_BOUNDED_SIDE,CGAL::ON_UNBOUNDED_SIDE>(other);
    
    if (queries.first == queries.second)
      throw  InvalidArgumentError("Surfaces must not collide.");
    else
    {
    
      return manipulate_vertex_selection_with_direction<0>(other,adjustment,max_iter);

    }   


}

/* -- Surface Mesh Properties -- */

/**
 * @brief Computes the average edge length in the stored mesh object. 
 * @param  none 
 * @return the average edge length
 */
inline double Surface::average_edge_length()
{
   double sum = 0; 
   for ( edge_descriptor e : mesh.edges())
   {
         halfedge_descriptor he = mesh.halfedge(e);
         sum+=static_cast<double>( CGAL::Polygon_mesh_processing::edge_length(he,mesh));
     
   }
   return sum/mesh.number_of_edges();
}

/* -- Surface Mesh Segmentation -- */ //TODO RENAME

/**
 * @brief Segments the surface, called Domain::boundary_segmentation function.
 *
 * Used CGAL function sharp edges segmentation that marks sharp edges 
 * and segments facet restricted by marked edges.   
 * @see [sharp_edges_segmentation](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__detect__features__grp.html) 
 *
 * @param nb_of_patch_plus_one used as the initial value to mark the surface segmentations. 
 * @return vector of facets points represented as triangles and associated tags. 
 */
inline  std::vector<std::pair<Surface::Triangle_3 , std::pair<int,int>>> Surface::surface_segmentation(int nb_of_patch_plus_one, double angle_in_degree)
{
    typedef boost::property_map<Mesh,CGAL::edge_is_feature_t>::type EIFMap; 
    EIFMap eif = get(CGAL::edge_is_feature, mesh);
    
    //std::size_t number_of_patches;                

    Mesh::Property_map<face_descriptor, std::pair<int,int> > patch_id_map;
    Mesh::Property_map<vertex_descriptor,std::set<std::pair<int,int> > > vertex_incident_patch_map;                      

    patch_id_map = mesh.add_property_map<face_descriptor,std::pair<int, int> >("f:pid",std::pair<int,int>()).first;
    vertex_incident_patch_map = mesh.add_property_map<vertex_descriptor,std::set<std::pair<int, int> > >("f:vip",std::set<std::pair<int, int> >()).first;
        
    CGAL::Polygon_mesh_processing::sharp_edges_segmentation(mesh, angle_in_degree, eif,patch_id_map,
                                     CGAL::Polygon_mesh_processing::parameters::first_index(nb_of_patch_plus_one)
                                    .vertex_incident_patches_map(vertex_incident_patch_map));

   std::vector<std::pair<Triangle_3,std::pair<int,int>>> Tri2tagvec;
   Vertex_point_pmap vpm = get(CGAL::vertex_point,mesh);

   for ( face_descriptor f : mesh.faces())
   {
      halfedge_descriptor he = mesh.halfedge(f);
      Point_3 p1 = get(vpm,mesh.source(he));
      he = mesh.next(he);
      Point_3 p2 = get(vpm,mesh.source(he));
      he = mesh.next(he);
      Point_3 p3 = get(vpm,mesh.source(he));

      Triangle_3 tri(p1,p2,p3);

      Tri2tagvec.push_back(std::make_pair(tri,get(patch_id_map,f)));     
   }
   return Tri2tagvec;
}

/**
 * @brief Checks if the surface bounds a volume.
 * @param none
 * @return true if surface bounds a volume  
 */
inline bool Surface::does_bound_a_volume()
{
   return CGAL::Polygon_mesh_processing::does_bound_a_volume(mesh);
}

/**
 * @brief Computes the distance between a point and the closest vertex in 
 * the triangulated surface.
 * 
 * @param point SVMTK Surface Point_3 object. 
 * @return the distance from point and the closest vertex in the surface
 */
inline double Surface::distance_to_point(Surface::Point_3 point)
{

   vertex_vector vertices = get_closest_vertices(point,1);
   if (vertices.size()==0)
       return 0;
   assert_non_empty_mesh();

   Vector_3 vector(mesh.point(vertices[0]),point);
  
   return CGAL::to_double(CGAL::sqrt(vector.squared_length()));
}

/**
 * @brief Handles the ouput for PreconditionError in boolean operations.
 *
 * @param other SVMTK Surface class object.
 * @return none 
 */
inline std::string Surface::CGAL_precondition_evaluation(Surface other)
{
      std::string output = "Following preconditions failed: "   ;       
      if (CGAL::Polygon_mesh_processing::does_bound_a_volume(mesh))
         output = output +"\n" +"- Surface does no bound a volume."  ;       
      if (CGAL::Polygon_mesh_processing::does_bound_a_volume(other.get_mesh()))
         output = output +"\n" + "- Argument surface does no bound a volume."  ;  
      if ( CGAL::Polygon_mesh_processing::does_self_intersect(mesh)) 
         output = output +"\n" + "- Surface self intersection in witihn union volume."  ;      
      if ( CGAL::Polygon_mesh_processing::does_self_intersect(other.get_mesh())) 
         output = output +"\n" + "- Argument surface self intersection in witihn union volume.";
         
      return output;
}

/**
 * @brief Sets the face oritentiaton to outwards 
 * @param none 
 * @return none  
 */
inline void Surface::set_outward_face_orientation()
{
  assert_non_empty_mesh();
  if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
  {
    CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
  }
}

/* -- Boolean Operations -- */
// ADD Preconditions Collision detectetion
/**
 * @brief Computes the intersection between two triangulated surface mesh.
 * 
 * Computes and corefine the intersection of two triangulated surfaces mesh,
 * @see [corefine_and_compute_intersection](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html) 
 * The functions have precondition that both surfaces bounds a volume 
 * and that both surfaces does not have self-intersections where the surfaces 
 * intersect.
 *
 *
 * @param other SVMTK Surface object
 * @return success true if intersection computation is successful  
 */
inline bool Surface::surface_intersection(Surface other)
{
   assert_non_empty_mesh();
   other.assert_non_empty_mesh();

   try
   {
      return CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh, other.get_mesh(),mesh);
   }
   catch (const std::exception &exc)
   {
      CGAL_precondition_evaluation( other);
      std::string output = "CGAL precondition error\n"+ CGAL_precondition_evaluation(other);
      throw PreconditionError(output.c_str());

   }     
   
}  

/**
 * @brief Computes the difference between two triangulated surface mesh.
 *
 * Computes and corefine the difference of two triangulated surfaces mesh,
 * @see [corefine_and_compute_difference](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html) 
 * The functions have precondition that both surfaces bounds a volume 
 * and that both surfaces does not have self-intersections.
 *
 *
 * @param other SVMTK Surface object
 * @return success true if intersection computation is successful  
 */
inline bool Surface::surface_difference(Surface other)
{    
   assert_non_empty_mesh();
   other.assert_non_empty_mesh();

   try 
   {
      return CGAL::Polygon_mesh_processing::corefine_and_compute_difference(mesh , other.get_mesh(), mesh);
   }
   catch (const std::exception &exc)
   {
      CGAL_precondition_evaluation(other);
      std::string output = "CGAL precondition error\n"+ CGAL_precondition_evaluation(other);
      throw PreconditionError(output.c_str());
   }   
}

/**
 * @brief Computes the union between two triangulated surface mesh.
 *  
 * Computes and corefine the union of two triangulated surfaces mesh,
 * @see [corefine_and_compute_union](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__corefinement__grp.html) 
 * The functions have precondition that both surfaces bounds a volume 
 * and that both surfaces does not have self-intersections in union 
 * volume.
 *
 * @param other SVMTK Surface object
 * @return success true if intersection computation is successful  
 */
inline bool Surface::surface_union(Surface other)
{  
   assert_non_empty_mesh();
   other.assert_non_empty_mesh();
   try 
   { 
       return CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh, other.get_mesh(), mesh); // CGAL::Polygon_mesh_processing::parameters::edge_is_constrained_map(is_constrained_map));
   }
   catch (const std::exception &exc)
   {
      std::string output = "CGAL precondition error\n"+ CGAL_precondition_evaluation(other);
      throw PreconditionError(output.c_str());
   }  
}
  
/**
 * @breif Keeps the largest connected component of the stored mesh.
 *
 * Computes and preserve the largest connected component in the stored mesh
 * @see [keep_largest_connected_components](https://doc.cgal.org/latest/Polygon_mesh_processing/group__keep__connected__components__grp.html) 
 *
 * @param none.
 * @return number of connected components.
 */ 
inline int Surface::keep_largest_connected_component()
{
 return CGAL::Polygon_mesh_processing::keep_largest_connected_components(mesh,1);
}

/**
 * @brief Finds a specified number mesh vertices that are closest to a point outside the mesh. 
 *
 * @param p1 is an arbiraty point outside the surface mesh 
 * @param num the maximum number of vertices that are added to the result
 * @return a vertex_vector that holds sum of `values`, or 0.0 if `values` is empty.
 */
inline Surface::vertex_vector Surface::get_closest_vertices(Point_3 p1, int num)
{
   
   Vertex_point_pmap vppmap = get(CGAL::vertex_point,mesh); 
   vertex_vector results;

   Tree tree(vertices(mesh).begin(),
            vertices(mesh).end(),
            Splitter(),
            Traits(vppmap)
   );

   Distance tr_dist(vppmap);
   
   K_neighbor_search search(tree,p1, num,0,true,tr_dist); 
 
   for ( auto vit=search.begin() ; vit!=search.end(); ++ vit) 
   {
        results.push_back(vit->first) ; 
   }
   return results;
} 

/**
 * @brief Finds a specified number mesh points that are closest to a point not on the mesh. 
 *
 * @param p1 is an arbiraty point not on the surface mesh 
 * @param num the maximum number of points that are added to the result
 * @return a point_vector that holds sum of `values`, or 0.0 if `values` is empty.
 */
inline Surface::point_vector Surface::get_closest_points(Point_3 p1, int num)
{
   vertex_vector vertices = get_closest_vertices( p1, num);
   return get_points(vertices);
} 

/** FIXME floating point error.
 * @brief Constructs a cylindric extension, which is a cylinder surface mesh combined with a sphere 
 * on the end closest to the mesh.
 * 
 * The user can use the union operation to combine it with the main mesh. 
 * The cylinder surface mesh is determined by centeriod of vertex points that are closest to a point 
 * outisde the mesh, radius, length and the option of crating a cylinder normal to the surface mesh. 
 * Works best on convex surfaces. 
 *
 * @param p1 is an arbiraty point outside the surface mesh, the closest vertex is 
 * @param radius of the cylinder surface
 * @param length of the cylinder surface from the surface mesh.
 * @param edge_length the edge length of the output surface.
 * @param normal if true the cylinder is set normal to the surface mesh, else in the direction of p1 
 * @return a SVMTK Surface class obejct.
 * 
 * @throws InvalidArgumentError if point is inside surface.
 */
inline std::shared_ptr<Surface> Surface::cylindric_extension(const Point_3& p1, double radius, double length, double edge_length, bool use_normal)
{
   assert_non_empty_mesh();
   
   if (  is_point_inside(p1) )
           throw  InvalidArgumentError("Select a point outside surface.");
   if (  length<=0 )
           throw  InvalidArgumentError("The length must be larger than 0.");         
           
   // FIXME 
   Point_3 p3, p4;
   FT normal_distance, point_distance; 
   Vector_3 normal_dir,point_dir;
   
   vertex_vector closest_vertex = get_closest_vertices(p1,1);
   vertex_vector vertices = get_closest_vertices(mesh.point(closest_vertex[0]),6); // used to be 8 

   std::vector<Point_3> points;
   for (auto vit : vertices)
   {
       points.push_back(mesh.point(vit));
   }

   Point_3 p2 = CGAL::centroid(points.begin(),points.end(),CGAL::Dimension_tag<0>()); 


   if ( p1 == p2 )  
     point_dir = CGAL::Polygon_mesh_processing::compute_vertex_normal(closest_vertex[0],mesh);
   else 
   {
     point_dir = p1-p2; 
     point_distance = CGAL::sqrt(point_dir.squared_length());
     if (point_distance==0)
     {
       std::cout << "#1" << std::endl;
     
     }
     point_dir=point_dir/point_distance;  
   }
   

   if (use_normal) 
   {    
        for (auto i : vertices ) 
        {
            normal_dir=normal_dir + CGAL::Polygon_mesh_processing::compute_vertex_normal(i,mesh);
        }
        normal_distance = CGAL::sqrt(normal_dir.squared_length());
        
        if (normal_distance==0)
          point_dir = point_dir;
        else  
          point_dir = normal_dir/normal_distance;  
   }
  
   p3 = p2 + point_dir*length;  
   p4 = p2 - point_dir*radius;
  
  /* double x0  = CGAL::to_double(p4.x());
   double y0  = CGAL::to_double(p4.y());
   double z0  = CGAL::to_double(p4.z());
   double x1  = CGAL::to_double(p3.x());
   double y1  = CGAL::to_double(p3.y());
   double z1  = CGAL::to_double(p3.z());
   */
   std::shared_ptr<Surface> cylinder(new Surface()); 

   cylinder->make_cylinder(p4,p3,radius,edge_length); 

   return cylinder; 
}

/**
 * @brief Constructs a cylindric connection bridge between the shortests line between two points 
 * in each surface.   
 *
 * The user can use the union operation to combine it with the main mesh. 
 * The cylinder surface mesh is determined by user 
 * arguments radius, length and the option of normal to the surface mesh or not. 
 * Works best on convex surfaces.
 *
 * @see Surface::assert_non_empty_mesh
 *
 * @param other SVMTK Surface class obejct.
 * @param radius of the cylinder surface
 * @param edge_length the edge length of the output surface.
 * @return a SVMTK Surface class object.
 *
 *
 */

inline std::shared_ptr<Surface> Surface::cylindric_connection(Surface other, double radius, double edge_length)
{
   assert_non_empty_mesh();
   Surface::vertex_vector results;
   Vertex_point_pmap vppmap = get(CGAL::vertex_point,other.get_mesh());  

   Tree tree(vertices(other.get_mesh()).begin(),
            vertices(other.get_mesh()).end(),
            Splitter(),
            Traits(vppmap)
   );
   Distance tr_dist(vppmap);
   FT distance;
   FT min_distance = FT(std::numeric_limits<double>::max());
   Point_3 query_point,p0,p1,p2,p3;
   Vector_3 normal;
   
   
   for (vertex_descriptor vit : mesh.vertices())
   {
        K_neighbor_search search(tree, mesh.point(vit), 2,0,true,tr_dist); 
        
        query_point = other.get_mesh().point((search.begin())->first);

        Point_3 current = mesh.point(vit);
        distance = CGAL::squared_distance(current,query_point);
        
        if (distance < min_distance)
        {
             min_distance = distance;
             p1 = query_point;
             p0 =  mesh.point(vit);
        }
   }
   normal = p1-p0; 
   FT length = CGAL::sqrt(normal.squared_length());
   normal=normal/length;  
  
   p2 = p0-normal*radius; 
   p3 = p1+normal*radius;
   /*
   double x0  = CGAL::to_double(p2.x());
   double y0  = CGAL::to_double(p2.y());
   double z0  = CGAL::to_double(p2.z());
   double x1  = CGAL::to_double(p3.x());
   double y1  = CGAL::to_double(p3.y());
   double z1  = CGAL::to_double(p3.z());
   */
   
   std::shared_ptr<Surface> cylinder(new Surface()); 
   cylinder->make_cylinder(p2,p3,radius,edge_length); 

   return cylinder;

}

/**
 * @brief Finds the shorest edge connected to a vertex for each vertex in the input.  
 *
 * Finds the shortest edge connectect to each vertex in the argument, and 
 * multiplies the edge with an adjustment factor. The adjustment is stored 
 * in a map with the vertices as a key.
 *
 * @param vector contains a subset of the surface mesh vertices. 
 * @param adjusmtent mulitplier of the shortest edge length. 
 * @return A map of the vertex and the movment of the vertex ( mostly used in relation of the vertex normal) 
 *
 * @throws Surface::assert_non_empty_mesh();
 */
inline std::map<Surface::vertex_descriptor, double> Surface::get_shortest_edge_map(const Surface::vertex_vector &vector, const double adjustment )
{

   assert_non_empty_mesh();
   CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
   
   std::map< vertex_descriptor, double> results;
   for (vertex_descriptor vit : vector)
   {
        Point_3 current = mesh.point(vit);
        HV_const_circulator vbegin(mesh.halfedge(vit),mesh), done(vbegin);
        FT min_edge = FT(100);
        do
        {
          FT temp = CGAL::sqrt(CGAL::squared_distance(current, mesh.point(*vbegin) )); 
          if (temp < min_edge )
          {
             min_edge = temp;
          }
        }while(++vbegin!=done);
       results[vit] = adjustment*static_cast<double>(CGAL::sqrt(min_edge));
   }
   return results;        
}

/**
 * @brief Finds the span in a given Cartesian direction, i.e. (min, max) vertex position
 * 
 * @param direction Cartesian direction ( 0=x ,1=y,2=z) 
 * @return a pair of doubles like (min,max)
 */
inline std::pair<double,double> Surface::span(int direction)
{
     assert_non_empty_mesh();
     auto bbox_3 = CGAL::Polygon_mesh_processing::bbox(mesh);
     std::pair<double,double> span (bbox_3.min(direction),bbox_3.max(direction) );
     return span;
}

/** 
 * @brief  Creates a surface mesh based on an implicit function
 * 
 * @tparam an implict function that takes Cartesian coordinates.
 *          The function has a boundary defined as  
 *         f(x,y,z)=0  and the interior defined as f(x,y,z) < 0 
 * 
 * @param implitict function  
 * @param bounding sphere that encloses the mesh construction. 
 * @param angular bound  
 * @param radius bound (TODO) 
 * @return void, mesh is constructed
 */
template<typename Implicit_function>
void Surface::implicit_surface(Implicit_function implicit_function,
             double bounding_sphere_radius,
             double angular_bound,
             double radius_bound,
             double distance_bound)
{
     surface_mesher(mesh,implicit_function,bounding_sphere_radius,angular_bound,radius_bound, distance_bound);
}

/** 
 * @brief Adjust the vertex coordinates of vertices in a vector that it iterated over.
 * 
 * @tparam Inputiterator 
 * @param begin call std::vector
 * @param end call of std::vector
 * @param mulitplier of the normal vector that determines the vertes movement.
 * @return void, mesh vertices are updated
 */
template<typename InputIterator >
void Surface::adjust_vertices_in_region(InputIterator begin , InputIterator  end, const double c)
{ 
  assert_non_empty_mesh();
  std::vector<std::pair<vertex_descriptor, Point_3> > adjust; 
  for ( ; begin != end; ++begin)
  {
      Vector_3 delta= CGAL::Polygon_mesh_processing::compute_vertex_normal(*begin,mesh); 
      Point_3 p = mesh.point(*begin) + c*delta;
      adjust.push_back(std::make_pair(*begin, p));
  }
  for (std::pair<vertex_descriptor, Point_3> s : adjust)
  {
      mesh.point(s.first) = s.second;
  }
}

/** TODO
 * @brief Adjust vertices in the map iterator input with the corresponding double value in the vertex normal direction.       
 *
 * @param map iterator begin  
 * @param map iterator end
 * @return void,  mesh is updated with new coordinates. 
 */
inline void Surface::adjust_vertices_in_region(vertex_scalar_map::iterator begin, vertex_scalar_map::iterator end) 
{
 assert_non_empty_mesh();
  std::vector<std::pair<vertex_descriptor, Point_3> > adjust; 
  for ( ; begin != end; ++begin)
  {
      Vector_3 delta= CGAL::Polygon_mesh_processing::compute_vertex_normal(begin->first,mesh); 
      Point_3 p = mesh.point(begin->first) + begin->second*delta;
      adjust.push_back(std::make_pair(begin->first, p));
  }
  for (std::pair<vertex_descriptor, Point_3> s : adjust)
  {
      mesh.point(s.first) = s.second;
  }
}

/** TODO
 * @brief Adjust vertices in the map iterator input with the corresponding vector direction.       
 *
 * @param map iterator begin  
 * @param map iterator end
 * @return void  class member mesh is updated with new coordinates for the vertices. 
 * @overload
 */
inline void Surface::adjust_vertices_in_region(vertex_vector_map::iterator begin, vertex_vector_map::iterator end) // map or two vectors
{
  for ( ; begin != end; ++begin)
  {
      Point_3 p = mesh.point(begin->first) + begin->second;
      mesh.point(begin->first) = p;
  }

}

/**
 * @brief Smooths the vertices in the iterator input. This is done by taking the sum of the vector edges for each vertex, and
 * multipled with the constant double input.  
 * 
 * 
 * @note removes all isolated vertices to avoid error in algorithm   
 * @tparam   
 * @param template iterator begin.
 * @param template iterator end.
 * @param const double c 
 * @return void; class member mesh is updated with new coordinates for the vertices. 
 */
template<typename InputIterator >
inline void Surface::smooth_laplacian_region(InputIterator  begin , InputIterator  end ,const double c)
{
  assert_non_empty_mesh();
  CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
  std::vector<std::pair<vertex_descriptor, Point_3> > smoothed;
  for ( ; begin != end; ++begin)
  {
      Point_3 current = mesh.point(*begin);
      Vector_3 delta=CGAL::NULL_VECTOR;
      HV_const_circulator vbegin(mesh.halfedge(*begin),mesh), done(vbegin);
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

/**
 * @brief Smooths the vertices in the iterator input. This is done by taking the sum of the vector edges for each vertex, and
 * multipled with the constant double input.  
 * 
 * @param
 * @param template iterator begin.
 * @param template iterator end.
 * @param const double c 
 * @return void; class member mesh is updated with new coordinates for the vertices. 
 */
inline void Surface::smooth_laplacian_region(vertex_vector_map::iterator begin, vertex_vector_map::iterator end ,const double c)
{
  assert_non_empty_mesh();
  CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
  std::vector<std::pair<vertex_descriptor, Point_3> > smoothed;
  for ( ; begin != end; ++begin)
  {
      Point_3 current = mesh.point(begin->first);
      Vector_3 delta=CGAL::NULL_VECTOR;
      HV_const_circulator vbegin(mesh.halfedge(begin->first),mesh), done(vbegin);
      do
      {
          delta += Vector_3(mesh.point(*vbegin) - current);
          *vbegin++;
      }while(vbegin!=done);

      Point_3 p = current + c*delta/mesh.degree(begin->first);  // WHY C
      smoothed.push_back(std::make_pair(begin->first, p));
  }
  for (std::pair<vertex_descriptor, Point_3> s : smoothed)
  {
      mesh.point(s.first) = s.second;
  }
}

/**
 * @brief Takes a CGAL polyhedron surface and copies the structure over to CGAL surface mesh.
 * 
 * 
 * @param polyhedron  defined CGAL::Polyhedron_3<Kernel> Polyhedron; 
 */
inline Surface::Surface(Polyhedron &polyhedron) 
{
    CGAL::copy_face_graph(polyhedron, mesh);
}

/**
 * @brief Read the surface corresponding to the input filename.
 * Current fileformats: 
           .off 
           .stl  
 * @param filename the string path to surface to load.
 */
inline Surface::Surface(const std::string filename)
{
    load_surface(filename,mesh);
}

/**
 * @brief Creates a SVMTK Surface class object usng points and connections 
 * between points.
 * @param points coordinates of vertices  
 * @param faces connections of  vertices
 */
inline Surface::Surface(std::vector<Point_3>& points,  std::vector<Face>& faces)
{
  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces);
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces,mesh);

  if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
  {
    CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
  }
}

/** 
 * @brief Takes a vector of verticrs and finds cluster of adjacent vertices with normals determined 
 * by  
 *       $ normal*n_2  > cos_angle $ 
 * with normal = Avg(n_i)   
 *
 * @note Smoothing causes most adjacent vertex normals to be consider within the cluster. Therefore, 
 * it is more robust to take the average of all vertex normal corresponding to the input.
 * Precondition of a non-empty mesh
 *  
 * @param[in,out] vertices a vector of vertices.  
 * @param cos_angle threshold cosinus angle between two vectors 
 * @return void 
 */
inline void Surface::get_normal_vector_cluster(Surface::vertex_vector &vertices, double angle_in_degree)  
{

  assert_non_empty_mesh();
  CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
  
  double cos_angle = std::cos(angle_in_degree*PI/180.0);

  int size = vertices.size();
  Vector_3 normal;
  
  // Compute the normal 
  for ( auto vit : vertices) 
  {
     normal = normal + CGAL::Polygon_mesh_processing::compute_vertex_normal(vit,mesh); 
  }
  normal = normal/CGAL::sqrt(normal.squared_length());
  
  for (int i = 0 ; i < size ; ++i)
  {
      HV_const_circulator vbegin(mesh.halfedge(vertices[i]),mesh), done(vbegin);
      do
      {
          Vector_3 n2 = CGAL::Polygon_mesh_processing::compute_vertex_normal(*vbegin,mesh);     
          
          if ( normal*n2 > abs(cos_angle))
          {
              vertex_descriptor tar =   *vbegin;
            
              if(std::find(vertices.begin(), vertices.end(), tar) == vertices.end()) 
              {
                 vertices.push_back(tar); 
                 ++size;
              } 
          } 
          *vbegin++;
      }while(vbegin!=done);
  }
}

/** 
 * @brief Taubin smothing of surface mesh. 
 *
 * Taubin smaoothing of the surface vertices. This corresponds to a 
 * Laplacian smoothing with value \lambda, followed by a Laplacian smoothing with value \mu
 * Given the requriemnt:     
 *               $ \lambda < -\mu $    
 * @note The Laplacian smoothin parameters are set, but the 
 *       user may consturct with their own parameters with smooth_laplacian function.
 * @param nb_iter number of iterations of smoothing   
 * @return void 
 */
inline void Surface::smooth_taubin(const size_t nb_iter) {
    for (size_t i = 0; i < nb_iter; ++i) {
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
 *       user combine smooth_laplacian.
 *
 * @param template iterator begin 
 * @param template iterator end 
 * @param nb_iter number of iterations of smoothing   
 * @return void 
 */
template<typename InputIterator >
void Surface::smooth_taubin_region(InputIterator begin , InputIterator end ,const size_t nb_iter)
{
    assert_non_empty_mesh();
    for (size_t i = 0; i < nb_iter; ++i) {
        this->smooth_laplacian_region(begin,end,0.8);
        this->smooth_laplacian_region(begin,end,-0.805);
    }
}

/** TODO : more functionallity
 * @brief Returns number of self intersection in the surface 
 * @param none    
 * @return number number of self intersection in the surface.
 */
inline int Surface::num_self_intersections() {
    assert_non_empty_mesh();
    std::vector< std::pair<face_descriptor, face_descriptor> > intersected_tris;
    CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(intersected_tris));
    return intersected_tris.size(); 
}

/**
 * @brief Combines smaller edges together, so that all edges are larger than the input parameter
 * @see [edge_collapse](https://doc.cgal.org/latest/Surface_mesh_simplification/index.html)
 * @param  target_edge_length the lower bound of the edges in the surface mesh.
 * @return r an interger that indicates the number of collapsed edges 
 */
inline int Surface::collapse_edges(const double target_edge_length)
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

/**
 * @brief Collapses smaller edges together, given the requirement that the edges have the same 
 * direction
 * @param none
 *   
 * @return r an interger that indicates the number of collapsed edges
 * @overload  
 */
inline int Surface::collapse_edges() {
    Cost_stop_predicate<Mesh> stop(1.e-6);
    const int r = CGAL::Surface_mesh_simplification::edge_collapse(
        mesh,
        stop);

    return r;
}

/**
 * @brief  Slices a surface mesh based on a plane, that is defined by the plane equation : 
 *       p1*x1 + p2*x2 +p3*x3 -x4 = 0.
 *
 * @note the ouput yields better result if the plane normal is close to unity. 
 * @param x1 plane equation parameter
 * @param x2 plane equation parameter
 * @param x3 plane equation parameter
 * @param x4 plane equation parameter
 * @return SVMTK slice object 
 */
template<typename Slice>
std::shared_ptr<Slice> Surface::mesh_slice(double x1,double x2, double x3 ,double x4)
{
   assert_non_empty_mesh();
   if ( x1==0 && x2==0 && x3==0)
      throw  InvalidArgumentError("Invalid plane parameters.");
   Plane_3 plane = Plane_3(x1, x2, x3, x4);
   return this->mesh_slice<Slice>(plane);
}

/**
 * @brief  Slices a SVMTK Surface class object according to a plane, 
 * and returns a SVMTK Slice class object.
 * TODO
 * @see [Plane_3] (https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html)   
 * @tparam SVMTK Slice class 
 * @param plane_3 a wrapped CGAL Plane_3 obejct 
 * @return slice SVMTK Slice class defined in Slice.h
 */
template<typename Slice>
std::shared_ptr<Slice> Surface::mesh_slice(Surface::Plane_3 plane_3) 
{
     assert_non_empty_mesh();
     
     typedef Kernel::Point_2 Point_2;
     typedef std::vector<Point_3>  Polyline_3; 
     typedef std::vector<Polyline_3> Polylines; 

     CGAL::Polygon_mesh_slicer<Mesh, Kernel> slicer(mesh); 
     Polylines polylines_3D;
     slicer(plane_3, std::back_inserter(polylines_3D));

     std::vector<std::vector<Point_2>> polylines_2;
     for ( auto  pol = polylines_3D.begin(); pol != polylines_3D.end(); ++pol ) 
     {
         std::vector<Point_2> result;
         std::vector<Point_2> polyline_2;
         for ( auto pit = pol->begin(); pit != pol->end(); ++pit)
         {
               polyline_2.push_back(plane_3.to_2d(*pit));
         }
         polylines_2.push_back(polyline_2);
     }
     std::shared_ptr<Slice> slice(new Slice(plane_3,polylines_2)); 
     return slice;
}  

/**
 * @brief  Slices a SVMTK Surface class object according to a plane, 
 * and returns the intersecting lines as a vetor of point vectors

 * @see [Plane_3] (https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html)   
 * @return slice SVMTK Slice class defined in Slice.h
 */
inline Surface::Polylines Surface::polylines_in_plane(Surface::Plane_3 plane_3) 
{
     assert_non_empty_mesh();
     CGAL::Polygon_mesh_slicer<Mesh, Kernel> slicer(mesh); 
     Polylines polylines;
     slicer(plane_3, std::back_inserter(polylines));

     return polylines;
}  


   

/**
 * @brief Isotropic remeshing of surface mesh. Remeshing of the surface mesh so 
 *        that all edges have the same length.
 *
 * 
 * CGAL isotropic_remeshing and split_long_edges.
 * @see [split_long_edges](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html) 
 * @see [isotropic_remeshing](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html) 
 * @note split_long_edges is used to avoid a pitfall described [here](https://doc.cgal.org/5.0.3/Polygon_mesh_processing/index.html)
 *
 * @param target_edge_length the edge length that is targeted in the remeshed patch. If 0 is passed then only the edge-flip, tangential relaxation, and projection steps will be done. 
 * @param nb_iter the number of iterations for the sequence of atomic operations performed. 
 * @param protect_border If true, constraint edges cannot be modified at all during the remeshing process. 
 * @return void 
 */
inline void Surface::isotropic_remeshing(double target_edge_length, unsigned int nb_iter, bool protect_border)
{
     assert_non_empty_mesh();
     CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
     CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh),
                              target_edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
                              .protect_constraints(protect_border));
}

/**
 * @brief Clips the surface mesh.
 * 
 * Clips a surface mesh based on a plane defined with the equation 
 *                    $ ax+by+cz +d = 0  $
 * @see [CGAL::Polygon_mesh_processing::clip](https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#Coref_section)
 * @see [CGAL::Plane_3](https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Plane__3.html)
 * @note Setting plane equation parameters close to unit values produce better clip.
 *
 * @param  a parameter in plane equation
 * @param  b parameter in plane equation
 * @param  c parameter in plane equation
 * @param  d parameter in plane equation
 * @param  preserve_manifold true to preserve manifold  
 * @return void 
 *
 */
inline bool Surface::clip(double a,double b, double c ,double d, bool preserve_manifold)
{      
   assert_non_empty_mesh();
   return CGAL::Polygon_mesh_processing::clip(mesh, Plane_3(a,b,c,d), CGAL::Polygon_mesh_processing::parameters::clip_volume(preserve_manifold));
}

/**
 * @brief Clips the surface mesh.
 * 
 * Clips a surface mesh given CGAL Plane_3 object defined by a point and a vector.
 * @param  point on the plane 
 * @param  vector a vector normal to the plane.
 * @param  preserve_manifold true to preserve manifold  
 * @return void 
 *
 * @overload
 */
inline bool Surface::clip(Point_3 point, Vector_3 vector , bool preserve_manifold)
{
   assert_non_empty_mesh();
   return CGAL::Polygon_mesh_processing::clip(mesh, Plane_3(point,vector), CGAL::Polygon_mesh_processing::parameters::clip_volume(preserve_manifold));
}

/**
 * @brief Clips the surface mesh.
 *
 * @param  point on the plane 
 * @param  vector a vector normal to the plane.
 * @param  radius of the circle 
 * @param  invert the clip
 * @param  preserve_manifold true to preserve manifold  
 * @return void 
 *
 * @overload
 */
inline bool Surface::clip(Point_3 point, Vector_3 vector, double radius, bool invert  , bool preserve_manifold)
{
   Surface circle;
   
   circle.make_circle_in_plane(point, vector,radius, radius/10. ) ; 
   return clip(circle,invert,preserve_manifold);
}

/**
 * @brief Clips a surface mesh given CGAL Plane_3 object.
 *
 * @param  plane 
 * @param  preserve_manifold true to preserve manifold  
 * @return void 
 * @overload 
 */
inline bool Surface::clip(Plane_3 plane, bool preserve_manifold)
{
   assert_non_empty_mesh();
   return CGAL::Polygon_mesh_processing::clip(mesh, plane, CGAL::Polygon_mesh_processing::parameters::clip_volume(preserve_manifold));
}

/**
 * @brief Clips a surface mesh given another SVMTK Surface object. 
 * @param  other SVMTK Surface object.
 * @param  preserve_manifold true to preserve manifold  
 * @param invert if true invert other surface.  
 * @return void 
 * @overload
 */
inline bool Surface::clip(Surface other,bool invert,bool preserve_manifold)
{
   assert_non_empty_mesh();
    if (invert) 
         CGAL::Polygon_mesh_processing::reverse_face_orientations(other.mesh);
    return CGAL::Polygon_mesh_processing::clip(mesh, other.mesh,CGAL::Polygon_mesh_processing::parameters::clip_volume(preserve_manifold));	

}

/**
 * @brief  Finds and fills holes in surface mesh.
 *
 * Uses CGAL function triangulate_refine_and_fair_hole. 
 * @see [triangulate_refine_and_fair_hole](https://doc.cgal.org/latest/Polygon_mesh_processing/group__hole__filling__grp.html)
 * Precondition non-empty mesh.
 * @param none.
 * @return nb_holes number of holes filled
 */
inline int Surface::fill_holes()
{
    //assert_non_empty_mesh();
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

/**
 * @brief  Triangulates faces of the surface mesh.
 *
 *
 * Uses CGAL function triangulate_faces.
 * @see [triangulate_faces] (https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html)
 * 
 * @param none    
 * @return true when done  
 */
inline bool Surface::triangulate_faces()
{
    assert_non_empty_mesh();
    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

    BOOST_FOREACH(face_descriptor fit, faces(mesh))
      if (next(next(halfedge(fit, mesh), mesh), mesh)
          !=   prev(halfedge(fit, mesh), mesh))
        std::cerr << "Error: non-triangular face left in mesh." << std::endl;
    return true;
}

/**
 * @brief Returns the points corresponding to vertices of the surfaces mesh.
 * @param vertices a vector of vertices, 
 * @return vector of surface points.
 */
inline std::vector<Surface::Point_3> Surface::get_points(Surface::vertex_vector &vertices) 
{
   assert_non_empty_mesh();
   point_vector result;
   for ( vertex_descriptor vit : vertices ){result.push_back(mesh.point(vit) );}
   return result;
} 

/**
 *  
 * @brief Returns the points of the surfaces mesh.
 * @param none   
 * @return vector of surface points.
 */
inline std::vector<Surface::Point_3>  Surface::get_points() 
{
   assert_non_empty_mesh();
   point_vector result;
   for ( vertex_descriptor vit : mesh.vertices()){result.push_back(mesh.point(vit) );}
   return result;
} 

/**
 * @brief Checks if a point is inside surface mesh.
 *  
 * Query if a point is inside the member surface mesh.
 * @param point_3 surface mesh point
 * @return bool true if point is inside surface otherwise false  
 */
inline bool Surface::is_point_inside(Point_3 point_3) 
{
     assert_non_empty_mesh();
     Surface::Inside is_inside_query(mesh);
     CGAL::Bounded_side res = is_inside_query(point_3);
     if (res == CGAL::ON_BOUNDED_SIDE or res == CGAL::ON_BOUNDARY)
         return true;
     else 
         return false;
}

/**
 * @brief Adjust surface mesh vertex coordinates in the normal vertex direction multiplied with an 
 * argument value.
 * 
 * @param c double value multiplier  
 * @return none
 */
inline void Surface::adjust_boundary(const double c)
{
    assert_non_empty_mesh();
    Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    Surface::adjust_vertices_in_region(vb, ve,c);
}

/** 
 * @brief Smooths all vertices of surface mesh using the smooth laplacian region function.  
 *
 * @see [Surface::smooth_laplacian_region]
 *
 * @param c a double multipler of the displacment vector that decides the new vertex coordinate.   
 * @param iter integer that decides the number of iterations of algorithm
 * @note zero and negative integer will have no effect, and should therfore be avoided.     
 * @return void, updated vertex coordinates 
 *  
 */
inline void Surface::smooth_laplacian(const double c, int iter)
{
    assert_non_empty_mesh();
    Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    for ( int i = 0 ; i< iter ; ++i)
    {
        Surface::smooth_laplacian_region(vb, ve,c);
    }
}

/**
 * @brief Smooth the triangulated surface. 
 * @see [CGAL::Polygon_mesh_processing::smooth_shape](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html)
 *
 * @param a time step that corresponds to the speed by which the surface is smoothed. A larger time step results in faster 
          convergence but details may be distorted to have a larger extent compared to more iterations with a smaller step. 
          Typical values scale in the interval (1e-6, 1].
 * @param nb_iterations number of iteations  
 * @return result vector of vertices  
 */
inline void Surface::smooth_shape(double time,int nb_iterations)
{
    assert_non_empty_mesh();
    CGAL::Polygon_mesh_processing::smooth_shape(mesh, time, CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iterations));
}

/** 
 * @brief Computes the centerline of the surface.
 * 
 * @see [CGAL::extract_mean_curvature_flow_skeleton](https://doc.cgal.org/latest/Surface_mesh_skeletonization/index.html)
 * 
 * @param none    
 * @return Polylines vector of points. 
 */
inline Surface::Polylines Surface::mean_curvature_flow() 
{
    assert_non_empty_mesh();
    
    typedef CGAL::Mean_curvature_flow_skeletonization<Surface::Mesh> Skeletonization;
    typedef Skeletonization::Skeleton                             Skeleton;
    typedef Skeleton::vertex_descriptor                           Skeleton_vertex;

    struct Construct_polylines
    {   
       const Skeleton& skeleton;
       std::vector<Surface::Point_3> polyline_3;
       std::vector<std::vector<Surface::Point_3>> polylines_3;
 
       Construct_polylines(const Skeleton& skeleton): skeleton(skeleton){}
       void start_new_polyline()
       {
           polyline_3.clear() ;
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
  
    std::vector<std::vector<Surface::Point_3>> polylines_3;
    Skeleton skeleton;
    CGAL::extract_mean_curvature_flow_skeleton(mesh, skeleton);
    Construct_polylines Visitor(skeleton);
    CGAL::split_graph_into_polylines(skeleton,Visitor);
    return Visitor.get_polylines();
}

/**
 * @brief Computes and returns the shortest surface path between two points. 
 *
 * Two points are projected on to the surface mesh, and the shortest surface path
 * is selected.
 *
 * @param x0 first point x coordinate
 * @param y0 first point y coordinate
 * @param z0 first point z coordinate
 * @param x1 second point x coordinate
 * @param y1 second point y coordinate
 * @param z1 second point z coordinate
 * @return a vector of points 
 */
inline Surface::Polyline Surface::get_shortest_surface_path(double x0, double y0, double z0, double x1, double y1, double z1)
{  
   Point_3 source(x0,y0,z0) ;
   Point_3 target(x1,y1,z1) ;
   return get_shortest_surface_path(source, target); 
}

/**
 * @brief Computes and returns the shortest surface path between two points.
 * @param source first point 
 * @param target second point.
 * @return a vector of points 
 * @overload 
 */
inline Surface::Polyline Surface::get_shortest_surface_path(Point_3 source, Point_3 target) 
{     
      assert_non_empty_mesh();

      typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Mesh> Traits;
      typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
      typedef Surface_mesh_shortest_path::Face_location Face_location;
      typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
      typedef CGAL::AABB_traits<Kernel, Primitive> AABB_Traits;
      typedef CGAL::AABB_tree<AABB_Traits> Tree;

      std::vector<Surface::Point_3> points;

      Tree tree(faces(mesh).first, faces(mesh).second, mesh); 
      Surface_mesh_shortest_path shortest_paths(mesh);     

      Face_location source_location = shortest_paths.locate(source,tree);
      Face_location target_location = shortest_paths.locate(target,tree);

      shortest_paths.add_source_point(source_location);
      shortest_paths.shortest_path_points_to_source_points(target_location.first,target_location.second , std::back_inserter(points));  
     
      return points;
} 

/**
 * @brief Creates a surface mesh structure with vertices and factes connecting vertices for a cube.
 * 
 * Constructes cube surface mesh with isotropic remeshin of each side.
 * 
 * @param x0 first point x coordinate
 * @param y0 first point y coordinate
 * @param z0 first point z coordinate
 * @param x1 second point x coordinate
 * @param y1 second point y coordinate
 * @param z1 second point z coordinate
 * @param edge_length the size of an edge in each direction
 * @return a vector of points 
 */
inline void Surface::make_cube( double x0, double y0, double  z0,  double x1, double y1, double z1, double edge_length) 
{
     typedef boost::multi_array<int, 3> array_type;
     
     std::vector<face_vector> sides;
     clear(); 
     if (x0==x1 or y0==y1 or z0==z1) 
        throw  InvalidArgumentError("Invalid argument.");



     double dx =(x1-x0);  
     double dy =(y1-y0);
     double dz =(z1-z0);
     
     
     if (dx<edge_length || dy<edge_length ||  dy<edge_length) 
        throw  InvalidArgumentError("Select smaller edge length.");    


     int Nx = static_cast<int>(dx/edge_length);
     int Ny = static_cast<int>(dy/edge_length);
     int Nz = static_cast<int>(dz/edge_length);
     
     double ddx = dx/static_cast<double>(Nx); 
     double ddy = dy/static_cast<double>(Ny); 
     double ddz = dz/static_cast<double>(Nz);  
      
     array_type map(boost::extents[Nx+1][Ny+1][Nz+1]);
         
     int index =0;
     for ( int i = 0 ; i< Nx+1 ; ++i)
     {
         for (int j=0 ; j < Ny+1 ; ++j )
         { 
              for (int k=0 ; k < Nz+1 ; ++k )
              {
                  if ( i==0 or i==Nx or j==0 or j==Ny or k==0 or k==Nz )
                  {     
                  mesh.add_vertex(
                       Point_3(x0+ static_cast<double>(i)*ddx,
                               y0+ static_cast<double>(j)*ddy,
                               z0+ static_cast<double>(k)*ddz) );
                     map[i][j][k] = index++;    
                  }     
              }
         }    
     }

     face_vector s1;
     face_vector s2;

     for (int j=0 ; j < Ny ; ++j )
     { 
         for (int k=0 ; k < Nz ; ++k )
         {
             s1.push_back(mesh.add_face(Index(map[0][j][k])  , Index(map[0][j+1][k+1]),Index(map[0][j+1][k])   ));
             s1.push_back(mesh.add_face(Index(map[0][j][k]), Index(map[0][j][k+1]) , Index(map[0][j+1][k+1])      ));
             s2.push_back(mesh.add_face(Index(map[Nx][j][k]) , Index(map[Nx][j+1][k]) , Index(map[Nx][j][k+1])  ));
             s2.push_back(mesh.add_face(Index(map[Nx][j+1][k]), Index(map[Nx][j+1][k+1]) , Index(map[Nx][j][k+1])   ));
         }
     }


     face_vector s3;
     face_vector s4;
     for (int i=0 ; i < Nx ; ++i )
     { 
         for (int k=0 ; k < Nz ; ++k )
         {
              s3.push_back(mesh.add_face(Index(map[i][0][k]) , Index(map[i+1][0][k])  ,Index( map[i][0][k+1])      ));
              s3.push_back(mesh.add_face(Index(map[i+1][0][k]),Index( map[i+1][0][k+1])  , Index(map[i][0][k+1])       ));
              s4.push_back(mesh.add_face(Index(map[i][Ny][k])  , Index(map[i][Ny][k+1]) , Index(map[i+1][Ny][k+1])   ));
              s4.push_back(mesh.add_face(Index(map[i+1][Ny][k]), Index(map[i][Ny][k]) , Index(map[i+1][Ny][k+1]) ));
         } 
     }

     face_vector s5;
     face_vector s6;
     for (int i=0 ; i < Nx ; ++i )
     { 
         for (int j=0 ; j < Ny ; ++j )
         {
              s5.push_back(mesh.add_face(Index(map[i][j][0] ), Index(map[i][j+1][0]) , Index(map[i+1][j+1][0]   ) ));
              s5.push_back(mesh.add_face(Index(map[i][j][0] ), Index(map[i+1][j+1][0])  , Index(map[i+1][j][0] )  ));
              s6.push_back(mesh.add_face(Index(map[i][j][Nz]), Index(map[i+1][j][Nz]), Index(map[i+1][j+1][Nz]  ) ));
              s6.push_back(mesh.add_face(Index(map[i][j][Nz]), Index(map[i+1][j+1][Nz]), Index(map[i][j+1][Nz])   ));
         } 
     } 
     
     sides.push_back(s1);
     sides.push_back(s2);
     sides.push_back(s3);
     sides.push_back(s4);
     sides.push_back(s5);
     sides.push_back(s6);     
     
     set_outward_face_orientation(); 

     for ( auto side : sides) 
     {     
     CGAL::Polygon_mesh_processing::isotropic_remeshing(side,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));        
     }
     return;
}

/**
 * @brief Creates a circle in a plane, with the plane defined by a point and vector
 *
 * @param point SVMTK Point_3 class object 
 * @param vector SVMTK Vector_3 class object 
 * @param radius of the circle
 * @param edge_length the target edge_length for the resulting circle.
 *
 *
 */
inline void Surface::make_circle_in_plane(Point_3 point, Vector_3 vector, double radius, double edge_length) 
{
     face_vector fv1,fv3;

     Index v0 = mesh.add_vertex(point);
     Index vb, vt;
   
     Plane_3 plane(point,vector);  

     Vector_3 t1 = plane.base1();
     Vector_3 t2 = plane.base2();
     
     t1=t1/std::sqrt(t1.squared_length()); 
     t2=t2/std::sqrt(t2.squared_length()); 
     
     
     int number_of_segments = static_cast<int>( 6.28*radius/edge_length);  
           
     if (number_of_segments<3) 
        throw  InvalidArgumentError("Select smaller edge length.");      
                   
     
     double c  = 360.0/(double)number_of_segments;
     
     for(int i = 0; i < number_of_segments; ++i) 
     {              
         Point_3 pb= point + t1*radius*std::cos(c*i*CGAL_PI/180) + t2*radius*std::sin(c*i*CGAL_PI/180); 
         vb = mesh.add_vertex(pb);

         if ( i!=0) 
             fv1.push_back(mesh.add_face(v0,vb,Index(vb-1)));

    }
    fv1.push_back(mesh.add_face(Index(0), Index(1), vb) );
    
    CGAL::Polygon_mesh_processing::isotropic_remeshing(fv1,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));

}     

/**
 * @brief Creates a surface mesh structure with vertices and factes connecting vertices for a cylinder
 * @param x0 x coordinate of the first cylinder center 
 * @param y0 y coordinate 0f the first cylinder center
 * @param z0 z coordinate 0f the first cylinder center
 * @param x1 x coordinate of the second cylinder center 
 * @param y1 y coordinate 0f the second cylinder center
 * @param z1 z coordinate 0f the second cylinder center
 * @param r0 the radius of the cylinder 
 * @param edge_length the target length of an edge in the mesh.
 * @return none
 * @overlaod
 */
inline void Surface::make_cylinder( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0,  double edge_length)
{    
     clear();
     Surface::make_cone( x0, y0,z0, x1,y1, z1, r0 , r0, edge_length ) ;
}

/**
 * @brief Creates a surface mesh structure with vertices and factes connecting vertices for a cone with
 * radius equal zero.
 *
 *
 * 
 * @param x0 x coordinate of the first cylinder center 
 * @param y0 y coordinate 0f the first cylinder center
 * @param z0 z coordinate 0f the first  cylinder center
 * @param x1 x coordinate of the second cylinder center 
 * @param y1 y coordinate 0f the second cylinder center
 * @param z1 z coordinate 0f the second cylinder center
 * @param r0 the none-zero radius of the cone, corresponding to the first cylinder center.  
 * @param edge_length the target length of an edge in the mesh.
 * @return none
 */
inline void Surface::make_cone_( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0, double edge_length) 
{
     clear();
     face_vector fv1,fv3;

     Index v0 = mesh.add_vertex(Point_3(x0,y0,z0));
     Index v1 = mesh.add_vertex(Point_3(x1,y1,z1));        
     Index vb, vt;  
  
     Vector_3 dir(Point_3(x0,y0,z0),Point_3(x1,y1,z1));
     double l1 = std::sqrt(dir.squared_length());
     Plane_3 plane(Point_3(x1,y1,z1),dir/l1   ); 
    
     Vector_3 t1 = plane.base1();

     Vector_3 t2 = plane.base2();
    
     t1=t1/std::sqrt(t1.squared_length()); 
     t2=t2/std::sqrt(t2.squared_length()); 
     
     int number_of_segments = static_cast<int>( 6.28*r0/edge_length);    
    
     if (number_of_segments<3) 
        throw InvalidArgumentError("Select smaller edge length.");    
        
     double c  = 360.0/(double)number_of_segments;
     
     for(int i = 0; i < number_of_segments; ++i) 
     {              
         Point_3 pb= Point_3(x0,y0,z0)  + t1*r0*std::cos(c*i*CGAL_PI/180) + t2*r0*std::sin(c*i*CGAL_PI/180); 
         vb = mesh.add_vertex(pb);

         if ( i!=0)
         {  
             fv1.push_back( mesh.add_face(v0,vb,Index(vb-1)) );
             fv3.push_back( mesh.add_face(Index(vb-1),vb,v1) );
         }

    }
    fv1.push_back(mesh.add_face(Index(2),vb,v0 ) );
    fv3.push_back(mesh.add_face(vb,Index(2),v1) );

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

/**
 * @brief Creates a surface mesh structure with vertices and factes connecting vertices for a cone.
 *
 * TODO
 *
 * @note The function also handles the special cases of:
 *                                                sharp cone  
 *                                                cylinder. 
 * @param x0 x coordinate of the first cone center 
 * @param y0 y coordinate 0f the first cone  center
 * @param z0 z coordinate 0f the first cone center
 * @param x1 x coordinate of the second cone center 
 * @param y1 y coordinate 0f the second cone center
 * @param z1 z coordinate 0f the second cone center
 * @param r0 the radius corresponding to the first cone center.   
 * @param r1 the radius corresponding to the second cone center. 
 * @param edge_length the target length of an edge in the mesh.
 * @return none 
 */
inline void Surface::make_cone( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0, double r1,  double edge_length) 
{


     if ((r0*r0 +r1*r1)==0) 
        throw  InvalidArgumentError("Select one non-zero radius.");   
         
     clear();
     if ( r0==0.0 )
     {
        make_cone_(x1,y1,z1,x0,y0,z0,r1,edge_length);
        return;
     }
     else if ( r1==0.0 )
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
     Plane_3 plane(Point_3(x1,y1,z1),dir/l1   );

     Vector_3 t1 = plane.base1();
     Vector_3 t2 = plane.base2();
     
     t1=t1/std::sqrt(t1.squared_length());    
     t2=t2/std::sqrt(t2.squared_length()); 
     
     double r = std::max(r0,r1);
     
     int number_of_segments = static_cast<int>( 6.28*r/edge_length); // approx pi, since setting it as int    
     
     if (number_of_segments<3) 
        throw  InvalidArgumentError("Select smaller edge length.");    
     
     double c  = 360.0/(double)number_of_segments;
     
     for(int i = 0; i < number_of_segments; ++i) 
     {         
         Point_3 pb= Point_3(x0,y0,z0)  + t1*r0*std::cos(c*i*CGAL_PI/180) + t2*r0*std::sin(c*i*CGAL_PI/180); 
         Point_3 pt = Point_3(x1,y1,z1)  + t1*r1*std::cos(c*i*CGAL_PI/180) + t2*r1*std::sin(c*i*CGAL_PI/180); 

         vb = mesh.add_vertex(pb);
         vt = mesh.add_vertex(pt);

         if ( i!=0)
         {  
             fv1.push_back( mesh.add_face(v0, vb,Index(vb-2)) );
             fv2.push_back( mesh.add_face(v1, Index(vt-2), vt) );
             fv3.push_back( mesh.add_face(Index(vt-2), Index(vb-2), vt) );
             fv3.push_back( mesh.add_face(vb, vt, Index(vb-2)) );
         }
    }
    fv1.push_back(mesh.add_face(Index(0), Index(2), vb) );
    fv2.push_back(mesh.add_face(Index(1), vt, Index(3) ) );

    fv3.push_back(mesh.add_face(vt, vb, Index(3))       );
    fv3.push_back(mesh.add_face(Index(2), Index(3), vb) );

                       

    CGAL::Polygon_mesh_processing::isotropic_remeshing(fv1,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));

    CGAL::Polygon_mesh_processing::isotropic_remeshing(fv2,
                              edge_length ,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
           
                              .protect_constraints(true));
                              
    CGAL::Polygon_mesh_processing::isotropic_remeshing(fv3,
                              edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(3)
                              .protect_constraints(true));
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
        static double function(double x, double y , double z) { return  (x-x0)*(x -x0) +  (y-y0)*(y -y0) + (z-z0)*(z -z0) -radius*radius; }
};


/**
 * @brief  FIXME ??
 * @note inline declaration of static variables in header files requrie c++ version equal or greate than 7 
 *       with the compiler option c++1z/c++17
 * @hideinitializer
 */
inline double sphere_wrapper::radius = 0; // inline declaration of static variable radius for sphere wrapper 
inline double sphere_wrapper::x0 = 0;     // inline declaration of static variable radius for sphere wrapper
inline double sphere_wrapper::y0 = 0;     // inline declaration of static variable radius for sphere wrapper
inline double sphere_wrapper::z0 = 0;     // inline declaration of static variable radius for sphere wrapper

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
 * @return void
 */
inline void Surface::make_sphere( double x0, double y0, double  z0,double r0, double edge_length) 
{    
   if (6.28*r0 < 3*edge_length) 
       throw  InvalidArgumentError("Select smaller edge length."); 
  
  
  sphere_wrapper sphere;
  sphere.radius = r0;
  sphere.x0=x0;
  sphere.y0=y0;
  sphere.z0=z0;
  surface_mesher(mesh,sphere.function,x0,y0,z0,r0,30,edge_length,edge_length);   

}

/**  
 * @brief splits edges to a target edge length
 * 
 * CGAL function for splitting long edges.
 * @see [split_long_edges](https://doc.cgal.org/latest/Polygon_mesh_processing/group__PMP__meshing__grp.html)
 * @param target_edge_length the maximum edge legnth that the longer edges can split into.
 * @return void 
 */
inline void Surface::split_edges(double  target_edge_length)
{
    assert_non_empty_mesh();
    
    CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
}

/**
 * @brief Saves the surface mesh to file.
 *
 * Valid file formats: off and stl. 
 * @param outpath string path to save file.
 *        @extensions : stl and off.
 * @return void
 */
inline void Surface::save(const std::string outpath)
{    
     assert_non_empty_mesh();
    
     std::string extension = outpath.substr(outpath.find_last_of(".")+1);
     std::ofstream out(outpath);
     if ( extension=="off")
     {
        out << mesh;
        out.close();
     } 
     else if ( extension=="stl")
     {
       write_STL(mesh,out);
     }
}

/**
 * @brief Reconstruct surface mesh. 
 * 
 * Reconstruct a surface based on a CGAL surface mesh object with points using CGAL poisson_reconstruction algorithm.
 * @see [poisson_reconstruction](https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/index.html)
 * @para bounding_sphere_radius indicating the radius of a sphere bounding the meshing operations
 * @para angular_bound bounds for the minimum facet angle in degrees.
 * @para radius_bound bound for the minimum for the radius of the surface Delaunay balls and the center-center distances respectively.
 * @para distance_bound bound for the minimum center-center distances respectively.
 * @overload  
 */
inline void Surface::reconstruct( double angular_bound, double radius_bound, double distance_bound )
{ 
    assert_non_empty_mesh();

    poisson_reconstruction(*this,angular_bound, radius_bound, distance_bound);
}






/**
 * @brief Computes the convex hull of the points.
 * @see [CGAL::convex_hull_3](https://doc.cgal.org/latest/Convex_hull_3/index.html)
 * @param none
 * @param result SVMTK Surface class.  
 */
inline std::shared_ptr< Surface > Surface::convex_hull()
{
    assert_non_empty_mesh();
    
    Polyhedron poly;
    auto point_vector = get_points(); 

    CGAL::convex_hull_3(point_vector.begin(), point_vector.end(), poly);
    auto result = std::make_shared< Surface >(Surface(poly));
    return result;
}

/** 
 * @brief  Returns a vector that contains the point and the corresponding computed vertex normal.
 * @param none 
 * @return a vector that contains the point and the corresponding computed vertex normal.
 */
inline std::vector<std::pair<Surface::Point_3, Surface::Vector_3>> Surface::get_points_with_normal() 
{
   assert_non_empty_mesh();
   
   std::vector<std::pair<Point_3, Vector_3>> result;
   
   for ( vertex_descriptor vit : mesh.vertices())
   {
        Vector_3 normal =  CGAL::Polygon_mesh_processing::compute_vertex_normal(vit,mesh);
        Point_3 point = mesh.point(vit); 
        result.push_back(std::make_pair(point,normal)  );
   }
   return result;
} 





#endif
