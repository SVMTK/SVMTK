// Copyright (C) 2018-2019 Lars Magnus Valnes
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
#ifndef __Slice_H


#define __Slice_H

/* --- Includes -- */
#include "SubdomainMap.h" 

/* -- STL -- */
#include <algorithm> 
#include <iterator>

/* -- CGAL 2D and 3D Linear Geometry Kernel -- */
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

/* -- CGAL 2D Conforming Triangulations and Meshes -- */
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

/* -- CGAL 2D Triangulation -- */
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>

/* -- CGAL Polyline_simplification_2 -- */
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/simplify.h> 
#include <CGAL/Polyline_simplification_2/Stop_below_count_ratio_threshold.h>

/* -- CGAL Bounding Volumes -- */
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_2.h>

/* -- CGAL IO -- */
#include <CGAL/IO/write_VTU.h>
#include <CGAL/IO/Triangulation_off_ostream_2.h>

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#define PI 3.14159265
//#include <CGAL/Arr_segment_traits_2.h>
//#include <CGAL/intersections.h>
/**
 * @brief Computes the length of an polyline 
 * i.e. vector of points
 * @tparam InputIterator iterator for a vector of points.
 * @param begin iterator for a vector of points 
 * @param end iterator for a vector of points   
 */
template< typename InputIterator> 
double length_polyline( InputIterator begin , InputIterator end)
{
  double length = 0.0;
  for(; begin != end; ++begin)
     length += static_cast<double>( CGAL::sqrt(CGAL::squared_distance(*begin, *(begin+1))));        
  return length; 
} 

/**
 * @brief Smooths a polygon, a closed loop of points, by minimizing the angles of the polygon.
 * @tparam Point_d, Vector_d 
 * @param polygon, a vector of points where the frist point is equal to the last. 
 * @param beta a smoothing factor.  
 * @return polyline which is smoothed
 */
template< typename Point_d, typename Vector_d> 
std::vector<Point_d> smooth_polygon( std::vector<Point_d>&  polygon, double beta)
{
    typename std::vector<Point_d> smoothed;
    Point_d ghost_p0 = polygon.rbegin()[1];
    Point_d ghost_pN = polygon.begin()[1];
    
    polygon.insert(polygon.begin(), ghost_p0); 
    polygon.push_back(ghost_pN);        
    for( typename std::vector<Point_d>::iterator pit = ++polygon.begin(); pit!=--polygon.end(); pit++) 
    {
         Point_d p0 =  *pit + beta*Vector_d(*pit,*(pit-1) ) + beta*Vector_d(*pit,*(pit+1));
         smoothed.push_back(p0);
    }
    return smoothed; 
}

/**
 * @brief Smooths a polyline by minimizing the angles of the polyline.
 * @tparam Point_d, Vector_d 
 * @param polyline  a vector of points 
 * @param beta a smoothing factor.  
 * @return polyline which is smoothed
 */
template< typename Point_d, typename Vector_d> 
std::vector<Point_d> smooth_polyline( std::vector<Point_d>& polyline, double beta)
{

    typename std::vector<Point_d> smoothed;
    
    if( polyline.size() < 3 ) 
        return polyline;    

    if( polyline.front() == polyline.back() ) 
        return smooth_polygon<Point_d,Vector_d>(polyline, beta); 

    smoothed.push_back(polyline.front());  
    for( typename std::vector<Point_d>::iterator pit = ++polyline.begin(); pit!=--polyline.end(); pit++) 
    {
         Point_d p0 =  *pit + beta*Vector_d(*pit,*(pit-1) ) + beta*Vector_d(*pit,*(pit+1));
         smoothed.push_back(p0);
    }
    smoothed.push_back(polyline.back() );
 
    return smoothed; 
} 

/**
 * @brief Checks if a 2D point is inside a polygon.
 * @tparam Point_2 
 * @param polygon a closed vector of points,
 * @param query 2D point 
 * @return true if point inside polygon, else it returns false. 
 */
template< typename Point_2> 
bool inside_polygon( std::vector<Point_2> polygon, Point_2 query ) 
{
     if( CGAL::ON_BOUNDED_SIDE == CGAL::bounded_side_2(polygon.begin(), polygon.end(), query )  )
         return true;
     else 
         return false; 
}

/**
 * @brief Computes the area of a point vector.
 * @tparam InputIterator iterator for a vector of points.
 * @param begin iterator for a vector of points 
 * @param end iterator for a vector of points   
 */
template< typename InputIterator> 
double area_of_facet_vector( InputIterator begin , InputIterator end)
{
  double area = 0.0;
  for(; begin != end; ++begin)
      area += static_cast<double>( CGAL::area(begin->vertex(0)->point(),
                                              begin->vertex(1)->point(),
                                              begin->vertex(2)->point()));        
  
  return area; 
}

/**
 * \struct
 * Used to store points and compute 
 * the minimum bounding radius required to 
 * enclose all of the added points 
 */
template< typename Kernel>
struct Minimum_sphere_2
{
    typedef typename CGAL::Min_sphere_of_spheres_d_traits_2< Kernel,typename Kernel::FT> Traits;
    typedef typename CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
    typedef typename Traits::Sphere                    Sphere;
    typedef typename std::vector<typename Kernel::Point_2> Polyline_2;
    typedef typename std::vector<Polyline_2> Polylines_2;  
  
    /**
     * @brief Adds vector of points to the struct.  
     * @param polyline a vector of points in 2D.
    */     
    void add_polyline( const Polyline_2 &polyline) 
    {
         for(auto it=polyline.begin();it != polyline.end(); ++it)
             S.push_back(Sphere(*it, 0.0));
    }
    
    /**
     * @brief Adds  vector of vector of points to the struct.  
     * @param polylines a vector of vector of points in 2D.
    */  
    void add_polylines(const Polylines_2 &polylines)
    {
         for(auto it=polylines.begin();it != polylines.end(); ++it)
         {
             for(auto pit=it->begin(); pit!=it->end(); ++pit)
                S.push_back(Sphere(*pit, 0.0));
         }
     } 
     
     /**
      * @brief Computes the minimum bounding radius required to enclose the added points.
      * @returns the radius that encloses all added points  
      */
     double get_bounding_sphere_radius()
     {
          Min_sphere ms(S.begin(), S.end());
           return CGAL::to_double(ms.radius());
     }
     private:
         std::vector<Sphere> S;
};

// DocString: Slice
/**
 * \Slice 
 * 
 * The SVMTK Slice class is used to store and manipulate triangulated surfaces in a plane, i.e. the thrid coordinate is neglected.  
 *
 * The SVMTK Slice class is implemnted with the Exact predicates inexact constructions kernel.
 * @see(CGAL::Exact_predicates_inexact_constructions_kernel)[https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html] 
 * 
 * @note The Slice does not handle cavities, but cavities can be assigned with adding surfaces and
 * optional SVMTK SubdomainMap object. @see Slice::add_surface_domains
 *              
 */
class Slice
{
    public :
       typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
 
       typedef Kernel::Plane_3          Plane_3;
       typedef Kernel::Point_2          Point_2;
       typedef Kernel::Point_3          Point_3;
       typedef Kernel::Vector_3         Vector_3;
       typedef Kernel::Segment_2        Segment_2;       
       typedef Kernel::Intersect_2      Intersect_2;
       typedef CGAL::Polygon_2<Kernel>  Polygon_2;
       typedef Kernel::FT FT; 
        // DelaunayMeshVertexBase_2` and `DelaunayMeshFaceBase_2`
       //typedef CGAL::Delaunay_mesh_vertex_base_2<Kernel>                Vbk;
       //typedef CGAL::Delaunay_mesh_face_base_2<Kernel>                  Fbk;
       typedef CGAL::Triangulation_vertex_base_with_info_2<int,Kernel>    Vb;
       //typedef CGAL::Triangulation_vertex_base_2<int,Kernel>    Vb;
       
       typedef CGAL::Delaunay_mesh_vertex_base_2<Kernel,Vb> Vbk;
              
       typedef CGAL::Triangulation_face_base_with_info_2<int,Kernel>      Fb_w_i;
       typedef CGAL::Constrained_triangulation_face_base_2<Kernel,Fb_w_i> C_fb_w_i;

       typedef CGAL::Delaunay_mesh_face_base_2<Kernel,C_fb_w_i> Fb;
    

       typedef CGAL::Triangulation_data_structure_2<Vbk, Fb> Tds;

       typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds,CGAL::Exact_predicates_tag> CDT1; //CGAL:No_constraint_intersection_requiring_constructions_tag
       
       typedef CGAL::Constrained_triangulation_plus_2<CDT1> CDT;
       
       typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
       typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher; // TODO test

       typedef CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold Stop;
       typedef CGAL::Polyline_simplification_2::Squared_distance_cost Cost;

       typedef CDT::Face_handle     Face_handle;
       typedef CDT::Vertex_handle   Vertex_handle;
       typedef CDT::Vertex_iterator Vertex_iterator;
       typedef CDT::Face_iterator   Face_iterator;
       typedef CDT::Edge            Edge;
       typedef CDT::Constraint_id   CID;
       
       typedef std::vector<Point_2>      Polyline_2;
       typedef std::vector<Polyline_2>   Polylines_2;   
       typedef std::vector<int>          Face;

      /** 
       *  stuct 
       *  @brief Used to sort polylines based on the number of points.
       */
       struct sort_vectors_by_size {
                  template<typename T>
                  bool operator()(const std::vector<T> & a, const std::vector<T> & b)
                       { return( a.size()>b.size() );}                   
       };
       
     /** 
      *  stuct 
      *  @brief Used to sort polylines based on the area FIXME 
      */
      struct sort_vectors_by_area {
                  template<typename T>
                  bool operator()(const std::vector<T> & a, const std::vector<T> & b)
                       { return ( area_of_facet_vector(a.begin(), a.end())>area_of_facet_vector(b.begin(), b.end()) );}                   
       };    

      /** 
      *  stuct 
      *  @brief Used to sort polylines based on the length of the point vector.
      */      
      struct sort_vectors_by_length {
                  template<typename T>
                  bool operator()(const std::vector<T> & a, const std::vector<T> & b)
                       { return ( length_polyline(a.begin(), a.end())>length_polyline(b.begin(), b.end()) );}                   
       };    
    
     // DocString: Slice     
     /** 
      * @brief Stores plane and polylines to member variables. 
      * @param plane 3D plane that represents the slice    
      */
      Slice(Plane_3 plane,Polylines_2 &polylines) : plane(plane)
      {
         min_sphere.add_polylines(polylines);
         constraints.insert( constraints.end(),polylines.begin(),polylines.end());
      }
 
      // DocString: Slice
      /** 
      * @brief Create an empty SVMTK Slice object.
      */
      Slice(){}
      
     // DocString: Slice
      /** 
      * @brief Constructs a SVMTK Slice object with a specified plane. 
      * @param plane plane of the slice. 
       
      */     
      Slice(Plane_3 plane) : plane(plane) {};
 
      // DocString: Slice
     /** 
      * @brief Constructs a SVMTK Slice object with a specified plane_3. 
      * @param point in the plane of the slice.  
      * @param vector normal to the plane of the slice.       
      */
      Slice(Point_3 point ,Vector_3 vector) : plane(point,vector) {};    
   
     // DocString: Slice  
     /** 
      * @brief Constructs a SVMTK Slice object with a specified plane_3. 
      *        The plane is defined with the equation 
      *                    $ ax+by+cz +d = 0  $
      * @param  a parameter in plane equation.
      * @param  b parameter in plane equation.
      * @param  c parameter in plane equation.
      * @param  d parameter in plane equation.
      */      
      Slice(double a, double b , double c , double d) : plane(a,b,c,d) {};

      ~Slice(){} 

     /** 
      * @brief Returns the minimum bounding radius to enclose all constraints.
      * @returns the minimum bounding radius to enclose all constraints.
      */
      double get_bounding_circle_radius()
      {    return min_sphere.get_bounding_sphere_radius();
      
      } 
      
      // DocString: set_plane      
      /** 
      * @brief Setz the plane member variable
      * @param inplane CGAL Plane_3 object.  
      */     
      void set_plane(Plane_3 plane)
      {   
         this->plane = plane;
      }
      
      // DocString: get_plane
      /** 
      * @brief Get the plane member variable
      * @returns the plane of SVMTK Slice object. 
      */      
      Plane_3& get_plane()
      {
         return this->plane;
      }
      
      // DocString: get_constraints
      /** 
      * @brief Get the constraints added to the class object
      * @returns the constraints added to the class object
      */     
      Polylines_2& get_constraints()
      {
          return constraints;
      }
     
     // DocString: clear_costraints
     /** 
      * @brief Clears all the constraints added to the class object
      */  
      void clear_costraints()
      { 
         constraints.clear();
      } 
     
     // DocString: number_of_subdomains
     /** 
      * @brief Returns the number of subdomains. 
      * @returns the number of subdomains
      */  
      int number_of_subdomains()
      {   
          return get_subdomains().size();
      } 
     
     // DocString: number_of_constraints
     /** 
      * @brief Returns the number of constraints 
      * @returns the number of constraints 
      */  
      int number_of_constraints()
      {   
         return constraints.size();
      }
      
     // DocString: number_of_faces
     /** 
      * @brief Returns the number of faces 
      * @return the number of faces, i.e triangles.
      */  
      int  number_of_faces()
      {   
         return cdt.number_of_faces();
      } 
     
     // DocString: get_edges
     /** 
      * @brief Returns the edges with the edge index.
      * @return map of edges and edge index.
      */  
      std::map<Edge,int>& get_edges()
      {  
          return this->edges;
      }    
      
     //DocString: get_points 
     /**
      * @brief Returns all points in the triangulation.
      * @return vector of points.
      */
      std::vector<Point_2> get_points()
      {
         std::vector<Point_2> out; 
         for(CDT::Point_iterator it = cdt.points_begin();  it!=cdt.points_end(); ++it)
            out.push_back(*it);
         
         return out;
      }

     //DocString: get_facet_points 
     /**
      * @brief Returns the mid points of each facet in the triangulation.
      * @return vector of points.
      */      
      std::vector<Point_2> get_facet_points()
      {
         std::vector<Point_2> out; 
         for(CDT::Face_iterator fit=cdt.faces_begin(); fit!=cdt.faces_end(); ++fit)
            out.push_back(CGAL::centroid( fit->vertex(0)->point(),fit->vertex(1)->point(),fit->vertex(2)->point() ));
            
         return out;
      }
      
     // DocString: remove_subdomain
     /** 
      * @brief Removes all faces in the mesh with specificed tags.  
      * @param tags vector of facet tags to be removed. 
      * @throws EmptyMeshError if cdt variable is empty.
      * @overlaod 
      */
      void remove_subdomain(std::vector<int> tags) 
      {
         assert_non_empty_mesh();
         for(CDT::Face_iterator fit=cdt.faces_begin(); fit!=cdt.faces_end(); ++fit)
         {
            if( std::find(tags.begin(), tags.end(), fit->info())!=tags.end() )
              cdt.delete_face(fit);
         } 
      }

      // DocString: remove_subdomain
     /**
      * @brief Removes all faces in the mesh with a specified tag   
      * @param tag removes facets with this integer tag
      * @overload
      */
      void remove_subdomain(int tag) 
      {
         std::vector<int> tags={tag}; 
         remove_subdomain(tags);
      }
      
      // DocString: create_mesh      
     /**  
      * @brief Create 2D mesh given a set of constraints, i.e. specifiec edges.
      * @param mesh_resolution a double value divisor used with the  
      *        minimum bounding radius to determine the maximum edge size.  
      */
      void create_mesh(double mesh_resolution) 
      {
          cdt.clear();
          double r = min_sphere.get_bounding_sphere_radius();
          double edge_size = r/mesh_resolution;
          std::cout<< "Edge size " << edge_size << std::endl;
          
          //remove_small_polygons( 10*edge_size ); -> FIXME 
          
          set_constraints(); 
          
          Mesher mesher(cdt);    
          mesher.set_criteria(Criteria(0.125, edge_size));

          std::cout << "Start  meshing."  << std::endl;
          mesher.refine_mesh(); 	
          //CGAL::lloyd_optimize_mesh_2(cdt,CGAL::parameters::number_of_iterations(10));
          std::cout << "Done  meshing."   << std::endl;
          
          
          std::map<Vertex_handle,bool> vertex_map;
          add_constraint_tags( ); 

          for(Vertex_iterator vit=cdt.vertices_begin(); vit!=cdt.vertices_end(); ++vit)       
          {
            if ( !cdt.is_infinite(vit) )
                 vertex_map[vit] = false;  
          }

          for(CDT::Face_iterator fit=cdt.faces_begin(); fit!=cdt.faces_end(); ++fit) // finite_faces_begin()
          {
             fit->info()=0;
             for(int i =0; i<3; ++i)
             {
                fit->vertex(i)->info()=0;
                if ( fit->is_in_domain() ) 
                    vertex_map[fit->vertex(i)] = true;
             }  
          }
          
          for(std::map<Vertex_handle, bool>::const_iterator it = vertex_map.begin();it != vertex_map.end(); ++it)   
          {

              if( !it->second) 
              { 
                if( cdt.are_there_incident_constraints(it->first) ) 
                    cdt.remove_incident_constraints(it->first); 
              }
          }
          cdt.remove_points_without_corresponding_vertex();
          remove_unbounded_structures();
      }
     
     /** 
      * @brief Create 2D mesh given a set of constraints, i.e. specifiec edges.
      * @param mesh_resolution a double value divisor used with the  
      *        minimum bounding radius to determine the maximum edge size.  
      */
      void create_mesh( double min_angle, double edge_size) 
      {
          cdt.clear();
          set_constraints(); 
 
          double B = 0.5/std::sin( min_angle*PI/180.0) ;

          Mesher mesher(cdt);    
          mesher.set_criteria(Criteria(B, edge_size));

          std::cout << "Start  meshing."  << std::endl;
          mesher.refine_mesh(); 	
          std::cout << "Done  meshing."   << std::endl;
          for(CDT::Face_iterator fit=cdt.faces_begin(); fit!=cdt.faces_end(); ++fit)
          {
             if( fit->is_in_domain() )  
                 fit->info()=0;
             else 
                 fit->info()=-1; // FIXME    
             for(int i =0; i<3; ++i)
                fit->vertex(i)->info()=0;
          }
      }

      // DocString: simplify
     /** 
      * @brief Simplify constraints to a specific point density.
      * @param point_density the number of points pr. average length unit 
      */
      void simplify(const double threshold=1.0)
      {       
          Polylines_2 result;
          std::cout<< "Start" << std::endl; 
          double accurate;
          for(auto c = this->constraints.begin(); c !=this->constraints.end(); ++c) 
          {   
             Polyline_2 temp; 
             accurate = std::max( threshold, 2./c->size() ) ;
             CGAL::Polyline_simplification_2::simplify(c->begin(), c->end(), Cost() , Stop(accurate), std::back_inserter(temp));
             result.push_back(temp); 
          }
          this->constraints.clear();
          this->constraints=result;
      }
 
     // DocString: connected_components
     /** 
      * @brief Calculates and returns the number of connected components.                                       
      * @returns number of connected components. 
      */
      int connected_components() 
      {
         assert_non_empty_mesh();
         std::vector<Face_handle> queue;
         std::map<Face_handle,bool> handled; 
         int num_cc=0;
         for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
            handled[fit] = false; 
         Face_handle fiq,fin;
         for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
         {
            if(handled[fit]) 
              continue;
            queue.push_back(fit);     
            while(!queue.empty())
            {
               fiq = queue.back(); 
               queue.pop_back();         
               if(handled[fiq]) 
                 continue;  
               handled[fiq] = true;
               for(int i =0; i < 3; i++) 
               {        
                  fin = fiq->neighbor(i); 
                  if(handled.find(fin) == handled.end()) 
                    continue;       
                  queue.push_back(fin);
               }
            }
            num_cc++;      
        }
        return num_cc;
     }

     // DocString: keep_largest_connected_component
     /** 
      * @breif Calculates and keeps the largest connected component.                                          
      */
      void keep_largest_connected_component() 
      {
          assert_non_empty_mesh();
          std::vector<Face_handle> queue;
          std::map<Face_handle,bool> handled; 
          std::vector<std::vector<Face_handle>> connected_components;
          int num_cc=0;
          for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
             handled[fit] = false; 
          Face_handle fiq,fin;
          for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
          {
             std::vector<Face_handle> connnected_component;
             if(handled[fit]) 
               continue;
             queue.push_back(fit);       
             while(!queue.empty())
             {
                fiq = queue.back(); 
                queue.pop_back();         
                if(handled[fiq]) 
                  continue;  
                handled[fiq] = true;
                connnected_component.push_back(fiq); 
                for(int i =0; i < 3; i++) 
                {        
                   fin = fiq->neighbor(i); 
                   if(handled.find(fin) == handled.end()) 
                     continue;     
                   queue.push_back(fin);
                }
             }
             num_cc++;
             connected_components.push_back(connnected_component); 
          }
          if( num_cc<2 ) 
            return;
          
          std::sort(connected_components.begin(), connected_components.end(), sort_vectors_by_size()); 
          for(auto ccit  = connected_components.begin()+1; ccit !=connected_components.end(); ccit++)
          {
             for(auto fit = ccit->begin(); fit!=ccit->end(); fit++)
                cdt.delete_face(*fit);
          } 
      }

      // DocString: add_constraints     
     /** 
      * @brief Adds constraints from another Slice object.
      * @param slice another SVMTK Slice object.
      */
      void add_constraints(Slice &slice) 
      {  
         add_constraints(slice.get_constraints());
      }
      
     // DocString: add_constraints      
     /** 
      * @brief Adds polylines to constraints.
      * @param polylines vector of vector of points in 2D.

      */
      void add_constraints(Polylines_2 &polylines) 
      {  
         for( auto polyline : polylines) 
              add_constraint(polyline); 
      }
           
     // DocString: add_connstraint
     /** 
      * @brief Adds polyline to constraints .
      * @param polyline vector of points in 2D.
      */
      void add_constraint(Polyline_2 &polyline) 
      {  
         std::vector<Point_2> update;
         Polylines_2 updates;
         if( !this->constraints.empty() )
         {  
            // Update constraints 
            for( auto line : this->constraints )
                updates.push_back(add_intersecting_point(line, polyline));       
                                
            // Update polyline
            Polyline_2 temp = polyline;  
            for( auto line : this->constraints )
                 temp = add_intersecting_point(temp, line); 

            updates.push_back(temp);      
                  
            min_sphere.add_polylines(updates);
            this->constraints.clear();
            this->constraints = updates;
 
         }            
         else 
         {
            min_sphere.add_polyline(polyline);
            this->constraints.push_back(polyline);
         }
 
      }
      
      /**
         TODO 
         @brief 
         
         @note
               
         @param area_threshold lower bound of polygon area.
         @return none                 
      
      */
      void remove_small_polygons(double area_threshold=10.)
      {
           for( auto cit = constraints.begin();  cit != constraints.end(); ) 
           {
               if( cit->back() == cit->front()  )
               {
                   auto polygon = Polygon_2(cit->begin(),cit->end());
                   auto area = abs(polygon.area());
                   if( area < FT(area_threshold)) 
                       cit = constraints.erase(cit); 
                   else 
                       cit++;                  
               }  
               else   
                  cit++; 
           }
      }
      

     /**
     TODO
       @brief 
       @param 
       @return 
     */
     std::vector<Point_2> add_intersecting_point(std::vector<Point_2> constraint, std::vector<Point_2> polyline )
     {
         std::vector<Point_2> update;
         for( auto cit1 = constraint.begin(), cit2 =  std::next(constraint.begin()); cit2 != constraint.end(); cit1++, cit2++)  
         {
            update.push_back(*cit1);
            for( auto pit1 = polyline.begin(), pit2 = std::next(polyline.begin()); pit2 != polyline.end(); pit1++, pit2++ )
            {
                auto res = find_intersection(*pit1,*pit2,*cit1,*cit2);
                if( res.first )
                {
                   if( update.front() != res.second ) 
                       update.push_back(res.second); 
                } 
            }
         }
         if( constraint.back() != update.back() ) 
            update.push_back( constraint.back() ) ;
         return update;
     }  
     
     /**
     TODO
       @brief 
       @param 
       @return 
     */
     std::pair<bool,Point_2> find_intersection(Point_2 p1,Point_2 p2, Point_2 q1 , Point_2 q2)
     {
         Segment_2 edge( q1, q2);
         Segment_2 query(p1, p2); 
         const auto result = CGAL::intersection(edge,query);
         if (result) 
         {
            const Point_2* p = boost::get<Point_2 >(&*result);
            return std::make_pair(true,Point_2(*p));
         }
         return std::make_pair(false, p1) ;
     } 
      
     /** 
      * @brief Adds constraints to the CGAL triangulation object cdt.
      * @param 
      * @return none 
      */
      void set_constraints() 
      {    
         for(auto pol: constraints)
         {
            if( !pol.empty() ) 
            {
               cids.push_back(cdt.insert_constraint(pol.begin(), pol.end()));
               this->bbox_2 += CGAL::bbox_2(pol.begin(), pol.end()); 
            }
         }
      }
      
     // DocString: slice_surfaces      
     /** 
      * @brief Slice a number of surfaces with the same plane, and store the constraints.
      * @tparam Surface SVMTK Surface object.
      * @param surfaces a vector of SVMTK Surface objects.
      */
      template<typename Surface> 
      void slice_surfaces(std::vector<Surface> surfaces) 
      {
         for(auto surf : surfaces) 
         {
            std::shared_ptr<Slice> temp = surf.template get_slice<Slice>(this->plane);              
            this->add_constraints(*temp.get()); 
         }
      }  
  
     // DocString: get_subdomains      
     /**
      * @brief Returns the subdomain tags in the mesh.   
      * @returns result a set of integers that represents the subdomain tags in the mesh.
      * @throws EmptyMeshError if cdt variable is empty.
      */
      std::set<int> get_subdomains()
      {
         assert_non_empty_mesh();
         std::set<int> result;
         for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
            result.insert(static_cast<int>(fit->info()));
         return result;
      }
   
     // DocString: export_as_surface      
     /** 
      * @brief Transforms the 2D mesh to 3D surface mesh stored in a SVMTK Surface object.
      * @precondition slice plane must be defined. 
      * @tparam Surface SVMTK Surface object.
      * @returns a SVMTK Surface object. 
      * @throws EmptyMeshError if cdt variable is empty.
      */
      template<typename Surface> 
      std::shared_ptr<Surface> export_as_surface() 
      {
         assert_non_empty_mesh(); 
         typedef std::vector<std::size_t> Face;
         typedef CDT::Vertex_handle Vertex_handle;

         std::vector<Point_3> points;
         std::vector<Face> faces;
         std::map<Vertex_handle, int> index_of_vertex;
         int i = 0;
         for(CDT::Point_iterator it = cdt.points_begin();  it!=cdt.points_end(); ++it, ++i)
         {
            points.push_back(plane.to_3d(*it));  
            index_of_vertex[it.base()] = i;
         }
         for(CDT::Face_iterator fit = cdt.faces_begin(); fit!=cdt.faces_end(); ++fit)
         {
            Face temp;
            temp.push_back( index_of_vertex[fit->vertex(0)]);
            temp.push_back( index_of_vertex[fit->vertex(1)]);
            temp.push_back( index_of_vertex[fit->vertex(2)]);
            faces.push_back(temp);
         }
         std::shared_ptr<Surface> surf(new Surface(points, faces)); 
         return surf;  
      }
      
     // DocString: export_as_surface      
     /** 
      * @brief Transforms the 2D mesh to 3D surface mesh stored in a SVMTK Surface object.
      * @precondition a vector of z coordinates
      * @tparam Surface SVMTK Surface object.
      * @returns a SVMTK Surface object. 
      * @throws EmptyMeshError if cdt variable is empty.
      */      
      template<typename Surface> 
      std::shared_ptr<Surface> export_as_surface( std::vector<double> z ) 
      {
         assert_non_empty_mesh(); 
         typedef std::vector<std::size_t> Face;
         typedef CDT::Vertex_handle Vertex_handle;

         std::vector<Point_3> points;
         std::vector<Face> faces;
         std::map<Vertex_handle, int> index_of_vertex;
         int i = 0;
         for(CDT::Point_iterator it = cdt.points_begin();  it!=cdt.points_end(); ++it, ++i)
         {
            points.emplace_back( it->x(), it->y() , z[i] );  
            index_of_vertex[it.base()] = i;
         }
         for(CDT::Face_iterator fit = cdt.faces_begin(); fit!=cdt.faces_end(); ++fit)
         {
            Face temp;
            temp.push_back( index_of_vertex[fit->vertex(0)]);
            temp.push_back( index_of_vertex[fit->vertex(1)]);
            temp.push_back( index_of_vertex[fit->vertex(2)]);
            faces.push_back(temp);
         }
         std::shared_ptr<Surface> surf(new Surface(points, faces)); 
         return surf;  
      }
  
     // DocString: save      
     /** 
      * @brief Saves the 2D mesh to file. 
      * Valid formats : off, stl, vtu and mesh(with tags).                                          
      * @param filename where the 2D mesh is to be stored.
      */
      void save(std::string outpath)
      {
         assert_non_empty_mesh(); 
         std::string extension = outpath.substr(outpath.find_last_of(".")+1);
         std::ofstream out(outpath);
         if( extension=="off" )
           CGAL::export_triangulation_2_to_off(out,cdt);
         else if( extension=="stl" )
           write_STL(outpath);
         else if( extension=="vtu" )
           CGAL::IO::write_VTU(out,cdt);  
         else if( extension=="mesh" )
           output_slice_to_medit_(out);
      }

     /** 
      * @brief Writes cdt to stl file.
      * @param filename where the 2D mesh is to be stored.
      * @return none 
      */
      void write_STL(const std::string filename)
      {
         std::ofstream file(filename);
         file.precision(6);
         file << "solid "<< filename << std::endl;
         for(CDT::Face_iterator fit=cdt.faces_begin(); fit!=cdt.faces_end(); ++fit)
         {
            Point_3 p1 = plane.to_3d( fit->vertex(0)->point());
            Point_3 p2 = plane.to_3d( fit->vertex(1)->point());
            Point_3 p3 = plane.to_3d( fit->vertex(2)->point());
       
            file << "facet normal " << plane.orthogonal_vector().x() << " " << plane.orthogonal_vector().x()    << " " << plane.orthogonal_vector().x() <<std::endl;
            file << "outer loop"<< std::endl;
            file << "\t" << "vertex " << p1.x() <<" "<< p1.y()   <<" "<< p1.z() << std::endl;  
            file << "\t" << "vertex " << p2.x() <<" "<< p2.y()   <<" "<< p2.z()   << std::endl; 
            file << "\t" << "vertex " << p3.x() <<" "<< p3.y()   <<" "<< p3.z()  << std::endl;             
            file <<"endloop" << std::endl;
            file <<"endfacet"<< std::endl;
         }
         file << "endsolid"  << std::endl;
      }

      /**
       * @brief Checks if the mesh is empty. 
       * @returns true if mesh is not empty. 
       * @throws EmptyMeshError if mesh is empty.
       */
      bool assert_non_empty_mesh()
      {  
         if( this->cdt.number_of_faces()==0 )
           throw EmptyMeshError("2D mesh object is empty.");
         return true;
      }
              
     /**
      * @brief Writes 2D mesh to medit file.
      *
      * Code based on CGAL output_to_medit, but for 2D meshes. Allows for edges to have a tag.
      * @see [output_to_medit] (https://github.com/CGAL/cgal/blob/master/Mesh_3/include/CGAL/IO/File_medit.h) 
      *
      * @param os ostream of the output file   
      * @throws EmptyMeshError if cdt variable is empty. 
      */
      void output_slice_to_medit_(std::ostream& os)
      {
         set_interface_tags(); 
         assert_non_empty_mesh(); 
         Tds tds = cdt.tds();
         os << std::setprecision(17);
         os << "MeshVersionFormatted 1\n"
            << "Dimension 2\n";

         os << "Vertices\n" << cdt.number_of_vertices() << '\n';
         boost::unordered_map<Vertex_handle, int> V;
         int inum = 1;
         for(Vertex_iterator vit=cdt.vertices_begin(); vit!=cdt.vertices_end(); ++vit)
         {
            V[vit] = inum++;
            vit->info() = 0; 
            os << *vit <<" "<< 0 <<'\n';
         }
         // FIXME 
         std::map<std::pair<Vertex_handle,Vertex_handle>,int> set_edges; 
         for(Face_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit) 
         {
            for(int i =0; i<3; ++i)
            {  
               Edge eit(fit,i);
               Vertex_handle vh1 = fit->vertex(cdt.ccw(i));
               Vertex_handle vh2 = fit->vertex(cdt.cw(i));
               if( V[vh1]>V[vh2] ) // Handles doubling of edges
                  set_edges[std::pair<Vertex_handle,Vertex_handle>(vh1,vh2)]=this->edges[eit];
               else
                  set_edges[std::pair<Vertex_handle,Vertex_handle>(vh2,vh1)]=this->edges[eit];  
            }
         }
         os << "Edges\n" 
            << set_edges.size() << '\n';
         for(auto eit : set_edges)      
            os << V[eit.first.first] << " " << V[eit.first.second]  <<" " << eit.second <<std::endl;         
         os << "Triangles\n" 
            << cdt.number_of_faces()<< '\n';
         for(Face_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit) 
            os << V[fit->vertex(0)]<< " " <<V[fit->vertex(1)]<<" " <<V[fit->vertex(2)]<< " " << fit->info() <<"\n";
         os << "End\n";
      }

     // DocString: add_surface_domains
     /** FIXME 
      * @brief Add tags to the facets in the 2D mesh based on overlapping surfaces and SubdomainMaps.  
      *                                    
      * Adds tags to faces dependent on position (inside/outside) according to closed 
      * triangulated surface in 3D.
      *  
      * @tparam Surface SVMTK Surface object.
      * @param surfaces a vector of SVMTK surface objects 
      * @param map derived from SVMTK AbstractMap objects, @see SubdomainMap.h 
      */
      template<typename Surface> 
      void add_surface_domains(std::vector<Surface> surfaces, AbstractMap& map) 
      {
         assert_non_empty_mesh();
         typedef boost::dynamic_bitset<>  Bmask;
         typedef std::pair<int,int>         Pid; 
         typedef std::map<Pid,int>      Pid_map;

         int fn,fi;
         Pid_map pid_map_;
         Pid spp;
         std::map<Face_handle,Bmask> masks;
         int index_counter=1;
         for(auto surf : surfaces) 
         {
            typename Surface::Inside inside(surf.get_mesh());
            for(CDT::Face_iterator fit=cdt.faces_begin(); fit!=cdt.faces_end(); ++fit)
            {
               Point_2 p2 = CGAL::centroid(fit->vertex(0)->point(),
                                           fit->vertex(1)->point(),
                                           fit->vertex(2)->point());  
               Point_3 p3 = this->plane.to_3d(p2);
               CGAL::Bounded_side res = inside(p3);
               if( res==CGAL::ON_BOUNDED_SIDE )
                   masks[fit].push_back(1);
               else
                   masks[fit].push_back(0);
            }
         } 
         for(auto bit : masks)
         {
            bit.first->info()=map.index(bit.second);
            if( bit.first->info()==0 )
              cdt.delete_face(bit.first);
         }
         Edge ei;

         for(Face_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit) 
         {
            fi = fit->info(); 
            for(int i =0;i<3;++i)
            {
               Edge ei(fit,i); 
               if( fit->neighbor(i)->is_in_domain() )
                 fn = fit->neighbor(i)->info();
               else 
                 fn = 0;    
              
               if( fn!=fi )
               {     
                 if( fn>fi )
                   spp={fn,fi};
                 else
                   spp={fi,fn};
                 std ::pair<Pid_map::iterator, bool> is_insert_successful = pid_map_.insert(std::make_pair(spp,index_counter));
                 if( is_insert_successful.second )
                   index_counter++;
                 this->edges[ei] = pid_map_[spp];    
               }
               else
                 this->edges[ei] = 0;        
            }
         }        
      } 
      
     // DocString: add_surface_domains
     /** FIXME 
      * @brief Add tags to the facets in the 2D mesh based on overlapping surfaces and SubdomainMaps.  
      *                                    
      * Adds tags to faces dependent on position (inside/outside) according to closed 
      * triangulated surface in 3D.
      *  
      * @tparam Surface SVMTK Surface object.
      * @param surfaces a vector of SVMTK surface objects 
      * @param map derived from SVMTK AbstractMap objects, @see SubdomainMap.h 
      */
      void set_interface_tags()
      {
           typedef std::pair<int,int>         Pid; 
           typedef std::map<Pid,int>      Pid_map;      
           if( !this->edges.empty() ) 
               return;
           int fn,fi;
           Pid_map pid_map_;
           Pid spp;      
           int index_counter=1;
           for(Face_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit) 
           {
               fi = fit->info();
               for(int i =0;i<3;++i)
               {
                  Edge ei(fit,i); 
                  if( fit->neighbor(i)->is_in_domain() )
                      fn = fit->neighbor(i)->info();
                  else 
                      fn = 0;    
              
                  if( fn!=fi )
                  {     
                      if( fn>fi )
                         spp={fn,fi};
                      else
                         spp={fi,fn};
                      std::pair<Pid_map::iterator, bool> is_insert_successful = pid_map_.insert(std::make_pair(spp,index_counter));
                      if( is_insert_successful.second )
                         index_counter++;
                      this->edges[ei] = pid_map_[spp];    
                  }
                  else
                      this->edges[ei] = 0;        
               }
           }             
      
      }

      void set_interface_tag(int subdomain, int tag)
      {
           if( !this->edges.empty() ) 
               return;
           
           for(Face_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit) 
           {
               if( fit->info()!=subdomain ) 
                  continue; 
               for(int i =0;i<3;++i)
               {
                  Edge ei(fit,i); 
                  if ( fit->neighbor(i)->info()!=subdomain)
                       this->edges[ei] = tag; 
               }
           }             
      
      }      

     /**
      * @brief Add triangles tags to triangles inside a given polygon.  
        @param polygon vector of 2D points with first and last element are the same 
        @param polygon_tag integer marker 
        @param replace_tag, overwriting a specific integer marker.  
      */      
      bool add_polygon_domain(std::vector<Point_2> polygon, int polygon_tag , int replace_tag=0) //TODO rename to overwrite? 
      {
           bool update = false;
           CGAL::Bbox_2 pbbox; 
           pbbox =  CGAL::bbox_2( polygon.begin() , polygon.end() );  
              
           if( !CGAL::do_overlap(this->bbox_2,pbbox))
               return false;
           
           update = true;
           for( CDT::Face_iterator fit=cdt.faces_begin(); fit!=cdt.faces_end(); ++fit )
           {  
               if( fit->info() != replace_tag) 
                   continue;
               Point_2 facet_point = CGAL::centroid( fit->vertex(0)->point(),
                                                     fit->vertex(1)->point(),
                                                     fit->vertex(2)->point() );
      
               if( CGAL::ON_BOUNDED_SIDE  == CGAL::bounded_side_2(polygon.begin(), polygon.end() , facet_point, Kernel()) ) 
                   fit->info() =  polygon_tag;
           }
           return update;
           
           
      }     
      
      
     // DocString: add_surface_domains       
     /**  
      * @brief Add tags to the facets in the 2D mesh based on overlapping surfaces and DefaultMap.                                      
      * @tparam Surface SVMTK Surface object.
      * @param surfaces a vector of SVMTK surface objects 
      * @overload 
      */
      template<typename Surface> 
      void add_surface_domains(std::vector<Surface> surfaces)
      {
           DefaultMap map =DefaultMap();
           add_surface_domains(surfaces,map);
      }  
      
     /**
      * @brief Returns vector of Faces.
      * @param None
      * @return vector of Faces  
      */
      std::vector<Face> get_facets() 
      {

         std::vector<Face> face_vector;
         boost::unordered_map<Vertex_handle, int> V;
         int inum = 1;
         
         for(Vertex_iterator vit=cdt.vertices_begin(); vit!=cdt.vertices_end(); ++vit)
            V[vit] = inum++;
         for(Face_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit) 
              face_vector.push_back( {V[fit->vertex(0)], V[fit->vertex(1)], V[fit->vertex(2)]} );
              
          return face_vector;
      }      
      
      /**
       * @brief Marks edges in each constraint with an unique tag  
       * @param None 
       * @return None 
       * 
       */
      void add_constraint_tags() 
      {
         int tag=0;
         Face_handle fit; 
 	 int i; 
         for( auto cid : this->cids ) 
         {    
              tag++;
              for( CDT::Vertices_in_constraint_iterator it1 = cdt.vertices_in_constraint_begin(cid), 
                   it2 = ++cdt.vertices_in_constraint_begin(cid); it2 != cdt.vertices_in_constraint_end(cid); ++it1,++it2 ) 
              {
                  if( cdt.is_edge(*it1,*it2, fit, i) ) // error ? 
	             this->edges[Edge(fit,i) ] = tag;
                  if( cdt.is_edge(*it2,*it1, fit, i) ) 
	             this->edges[Edge(fit,i) ] = tag;
              }
 	  }   
      }


      /**
       * @brief TODO FIXME
       * @param none 
       * @return vertex_indicies for all features.
       */
      std::vector<std::vector<int>> get_feature_vertices() 
      {
    
 	  std::vector<std::vector<int>> lines;
 	  boost::unordered_map<Vertex_handle, int> V;
         int inum = 0;
         for(Vertex_iterator vit=cdt.vertices_begin(); vit!=cdt.vertices_end(); ++vit)
         {
            V[vit] = inum++;
         }
 	 
         for( auto cid : this->cids ) 
         {    
              std::vector<int> line;
              
              for( CDT::Vertices_in_constraint_iterator vit = cdt.vertices_in_constraint_begin(cid);
                                                      vit != cdt.vertices_in_constraint_end(cid); ++vit) 
              {
                  if ( V.find(*vit) != V.end() )
                       line.push_back(V[*vit]);
              }    
              lines.push_back(line);
 	  }   
 	  return lines;
      }


      /**
       * @brief Removes triangles not bounded by constraints 
       * @param none 
       * @return none 
       */
      void remove_unbounded_structures() 
      {
          std::map<Vertex_handle,bool> vertex_map;
          // mark all vertices
          for(Vertex_iterator vit=cdt.vertices_begin(); vit!=cdt.vertices_end(); ++vit)       
          {
            if ( !cdt.is_infinite(vit) )
                 vertex_map[vit] = false;  
          }
          for(CDT::Face_iterator fit=cdt.faces_begin(); fit!=cdt.faces_end(); ++fit) // finite_faces_begin()
          {
             if( fit->is_in_domain() ) 
                for(int i =0; i<3; ++i)
                    vertex_map[fit->vertex(i)] = true;
             else 
                cdt.delete_face(fit);
          }
          for( auto vit : vertex_map ) 
          {
             if( !vit.second ) 
                cdt.delete_vertex(vit.first);
          }
      }
      
      
    private:
       boost::unordered_map<Face_handle,double>  triangle_data  ;
       boost::unordered_map<Vertex_handle,double>  point_data ;
       Minimum_sphere_2<Kernel> min_sphere;
       Polylines_2 constraints;
       CDT cdt;
       std::vector<CID> cids; //EXPERIMENTAL 
       Plane_3 plane;
       CGAL::Bbox_2 bbox_2; 
       std::map<Edge,int> edges; 

};

#endif 
