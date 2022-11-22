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

/* -- CGAL Polyline_simplification_2 -- */
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/simplify.h> 
#include <CGAL/Polyline_simplification_2/Stop_below_count_ratio_threshold.h>

/* -- CGAL Bounding Volumes -- */
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_2.h>

/* -- CGAL IO -- */
//#include <CGAL/IO/write_off_points.h>
//#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_VTU.h>
#include <CGAL/IO/Triangulation_off_ostream_2.h>

#define PI 3.14159265

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
 * @brief Computes the length of an polyline 
 * i.e. vector of points
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
 
       typedef Kernel::Plane_3   Plane_3;
       typedef Kernel::Point_2   Point_2;
       typedef Kernel::Point_3   Point_3;
       typedef Kernel::Vector_3   Vector_3;

       typedef CGAL::Triangulation_vertex_base_with_info_2<int,Kernel>    Vb;
       typedef CGAL::Triangulation_face_base_with_info_2<int,Kernel>      Fb_w_i;

       typedef CGAL::Constrained_triangulation_face_base_2<Kernel,Fb_w_i> C_fb_w_i;

       typedef CGAL::Delaunay_mesh_face_base_2<Kernel,C_fb_w_i> Fb;

       typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;

       typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds,CGAL::Exact_predicates_tag> CDT1; 
       typedef CGAL::Constrained_triangulation_plus_2<CDT1> CDT;
       
       typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
       typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

       typedef CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold Stop;
       typedef CGAL::Polyline_simplification_2::Squared_distance_cost Cost;

       typedef CDT::Face_handle     Face_handle;
       typedef CDT::Vertex_handle   Vertex_handle;
       typedef CDT::Vertex_iterator Vertex_iterator;
       typedef CDT::Face_iterator   Face_iterator;
       typedef CDT::Edge            Edge;

       typedef std::vector<Point_2> Polyline_2;
       typedef std::vector<Polyline_2> Polylines_2;   


       /** 
        *  stuct 
        *  @brief Used to sort polylines based on the number of points.
        */
       struct sort_vectors_by_size {
                  template<typename T>
                  bool operator()(const std::vector<T> & a, const std::vector<T> & b)
                       { return( a.size()>b.size() );}                   
       };

      struct sort_vectors_by_area {
                  template<typename T>
                  bool operator()(const std::vector<T> & a, const std::vector<T> & b)
                       { return ( area_of_facet_vector(a.begin(), a.end())>area_of_facet_vector(b.begin(), b.end()) );}                   
       };    
     

      
      struct sort_vectors_by_length {
                  template<typename T>
                  bool operator()(const std::vector<T> & a, const std::vector<T> & b)
                       { return ( length_polyline(a.begin(), a.end())>length_polyline(b.begin(), b.end()) );}                   
       };    
     
     
     
     
       
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
          set_constraints(); 
          double r = min_sphere.get_bounding_sphere_radius();
          double edge_size = r/mesh_resolution;
          std::cout<< "Edge size " << edge_size << std::endl;
          Mesher mesher(cdt);    
          mesher.set_criteria(Criteria(0.125, edge_size));

          std::cout << "Start  meshing."  << std::endl;
          mesher.refine_mesh(); 	
          std::cout << "Done  meshing."   << std::endl;
          for(CDT::Face_iterator fit=cdt.faces_begin(); fit!=cdt.faces_end(); ++fit)
          {
             fit->info()=0;
             for(int i =0; i<3; ++i)
                fit->vertex(i)->info()=0;
          }
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
             fit->info()=0;
             for(int i =0; i<3; ++i)
                fit->vertex(i)->info()=0;
          }
      }

      // DocString: simplify
     /**  
      * @brief Simplify constraints to a specific point density.
      * @param point_density the number of points pr. length
      */
      void simplify(const double point_density=0.4)
      {       
          Polylines_2 result;
          for(auto c = this->constraints.begin(); c !=this->constraints.end(); ++c) 
          {
             Polyline_2 temp; 
             double length = length_polyline(c->begin(),c->end()); 
             double adjustment = point_density*length/(double)(c->size());    
             CGAL::Polyline_simplification_2::simplify(c->begin(), c->end(), Cost() , Stop(adjustment), std::back_inserter(temp));
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
          // TODO: test sorting by area. 
          
           
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
         min_sphere.add_polylines(polylines);
         constraints.insert(constraints.end(), polylines.begin(), polylines.end());
      }
           
     // DocString: add_connstraint
     /** 
      * @brief Adds polyline to constraints .
      * @param polyline vector of points in 2D.
      */
      void add_constraint(Polyline_2 &polyline) 
      {  
         min_sphere.add_polyline(polyline);
         constraints.push_back(polyline);
      }
      

     /** 
      * @brief Adds constraints to the CGAL triangulation object cdt.
      */
      void set_constraints() 
      {    
         for(auto pol: constraints)
            cdt.insert_constraint(pol.begin(), pol.end());
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
      * Based on CGAL output_to_medit, but for 2D meshes. Allows for edges to have a tag.
      * @see [output_to_medit] (https://github.com/CGAL/cgal/blob/master/Mesh_3/include/CGAL/IO/File_medit.h) 
      *
      * @param os ostream of the output file   
      * @throws EmptyMeshError if cdt variable is empty. 
      */
      void output_slice_to_medit_(std::ostream& os)
      {
         assert_non_empty_mesh(); 
         Tds tds = cdt.tds();
         os << std::setprecision(17);
         os << "MeshVersionFormatted 1\n"
            << "Dimension 2\n";
         //-------------------------------------------------------
         // Vertices
         //-------------------------------------------------------
         os << "Vertices\n" << cdt.number_of_vertices() << '\n';
         boost::unordered_map<Vertex_handle, int> V;
         int inum = 1;
         for(Vertex_iterator vit=cdt.vertices_begin(); vit!=cdt.vertices_end(); ++vit)
         {
            V[vit] = inum++;
            vit->info() = 0; 
            os << *vit <<" "<< 0 <<'\n';
         }
         std::map<std::pair<Vertex_handle,Vertex_handle>,int> set_edges; 
         for(Face_iterator fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit) 
         {
            for(int i =0; i<3; ++i)
            {  
               Edge eit(fit,i);
               // TODO: check fit->vertex 
               Vertex_handle vh1 = fit->vertex(cdt.ccw(i));
               Vertex_handle vh2 = fit->vertex(cdt.cw(i));
               if( V[vh1]>V[vh2] )
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
            os << V[fit->vertex(0)]<< " " <<V[fit->vertex(1)]<<" " <<V[fit->vertex(2)]<< " " <<fit->info() <<"\n";
         os << "End\n";
      }

     // DocString: add_surface_domains
     /**
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
         typedef boost::dynamic_bitset<>   Bmask;
         typedef std::pair<int,int>  Pid; 
         typedef std::map<Pid,int> Pid_map;
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
               Point_2 p2 = CGAL::centroid(fit->vertex(0)->point() ,fit->vertex(1)->point(),fit->vertex(2)->point());  
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
    private:
       Minimum_sphere_2<Kernel> min_sphere;
       Polylines_2 constraints;
       CDT cdt;
       Plane_3 plane;
       std::map<Edge,int> edges; 

};

#endif 
