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
#include "SubdomainMap.h" 

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_conformer_2.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/simplify.h> 
#include <CGAL/Polyline_simplification_2/Stop_below_count_ratio_threshold.h>

#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/centroid.h>

#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_vtu.h>
#include <CGAL/IO/Triangulation_off_ostream_2.h>

#include <CGAL/utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/assertions.h>

#include <boost/unordered_map.hpp>

#include <algorithm> 
#include <iterator>
#include <queue>
#include <assert.h>  
#include <iterator>
#include <fstream>

/**
 * Function to compute the length of an polyline 
 * i.e. vector of points
 * @param begin iterator for an vector of points 
 * @param end iterator for an vector of points   
 */
template< typename InputIterator> 
double length_polyline( InputIterator begin , InputIterator end)
{
  double length = 0.0;
  for ( ; begin != end; ++begin)
  {
     length += static_cast<double>( CGAL::sqrt(CGAL::squared_distance(*begin, *(begin+1) ) ) ); 
  }         
  return length; 
} 

/**
 * @struct:
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
     * Adds vector of points to the struct  
     * @param polyline a vector of points in 2D.
     * @return none 
    */     
    void add_polyline( const Polyline_2 &polyline) 
    {
         for (auto it=polyline.begin();it != polyline.end(); ++it)
         {
             S.push_back(Sphere(*it, 0.0));
         }
    }
    
    /**
     * Adds  vector of vector of points to the struct  
     * @param polylines a vector of vector of points in 2D.
     * @return none 
    */  
    void add_polylines(const Polylines_2 &polylines)
    {
         for (auto it=polylines.begin();it != polylines.end(); ++it)
         {
             for ( auto pit=it->begin(); pit!=it->end(); ++pit)
             { 
                 S.push_back(Sphere(*pit, 0.0));
             }
         }
     } 
     
     /**
      * Computes the minimum bounding radius required to enclose the added points.
      * @param none   
      * @return the radius that encloses all added points  
      */
     double get_bounding_sphere_radius()
     {
          Min_sphere ms(S.begin(), S.end());
           return CGAL::to_double(ms.radius());
     }
     private:
         std::vector<Sphere> S;
};

/**
 *  
 */
class Slice
{
    public :
       // See (https://doc.cgal.org/latest/Kernel_23/group__kernel__predef.html)
       typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
 
       typedef Kernel::Plane_3   Plane_3;
       typedef Kernel::Point_2   Point_2;
       typedef Kernel::Point_3   Point_3;

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

       Slice(Plane_3 plane ,Polylines_2 &polylines);
       Slice(Plane_3 plane_3) : plane(plane_3) {};
       Slice(){}
       ~Slice(){} 

       void set_plane(Plane_3 inplane){ this->plane = inplane;}
       Plane_3& get_plane(){return this->plane;}
       template<typename Surface> 
       void add_surface_domains(std::vector<Surface> surfaces, AbstractMap& map); //todo rename
       template<typename Surface> 
       void add_surface_domains(std::vector<Surface> surfaces); //todo rename
       template< typename Surface> 
       void slice_surfaces(std::vector<Surface> surfaces);

       void remove_subdomains(std::vector<int> tags); 
       void remove_subdomains(int tag);

       void create_mesh(double mesh_resolution);    
       void simplify( const double point_density=0.4 );
       int  connected_components(); 
       void keep_largest_connected_component(); 

       void add_constraints(Polylines_2 &polylines); 
       void add_constraints(Slice &slice);
       void add_constraint(Polyline_2 &polyline);

       Polylines_2& get_constraints() { return constraints;}

       void clear_costraints(){ constraints.clear(); return; } 
       void set_constraints();

       int  number_of_subdomains(){return get_subdomains().size();} 
       int  number_of_constraints() { return constraints.size();}
       int  number_of_faces(){return cdt.number_of_faces();} 

       std::set<int> get_subdomains();
       std::map<Edge,int>& get_edges(){return this->edges;}  

       template<typename Surface> 
       std::shared_ptr<Surface> as_surface();

       double get_bounding_circle_radius(){return min_sphere.get_bounding_sphere_radius();} 
       void save(std::string outpath);
       void output_slice_to_medit_(std::ostream& os);
       void write_STL(const std::string filename);

       struct sort_vectors_by_size {
                  template<typename T>
                  bool operator()(const std::vector<T> & a, const std::vector<T> & b)
                       { return a.size() > b.size(); }                   
       };
    private:
       Minimum_sphere_2<Kernel> min_sphere;
       Polylines_2 constraints;
       CDT cdt;
       Plane_3 plane;
       std::map<Edge,int> edges; 

};

/**
 * Returns a set of integer that represents the subdomains in the mesh.   
 * 
 * @param none 
 * @return result a set of integers that represents the subdomain tags in the mesh.
 */
inline std::set<int> Slice::get_subdomains()
{
   std::set<int> result;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
       result.insert( static_cast<int>(fit->info()));
   }
   return result;
}

/**
 * Removes all cells in the mesh with a specified integer tag   
 * @param tag removes cells with this integer tag
 * @overload
 */
inline void Slice::remove_subdomains(int tag) 
{
    std::vector<int> tags; 
    tags.push_back(tag);
    remove_subdomains(tags);
}

/**
 * Removes all cells in the mesh with tags in a vector.  
 * @param tags vector of cell tag to be removed. 
 * @overload
 */
inline void Slice::remove_subdomains(std::vector<int> tags) 
{
    for(CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
    {
        if(std::find(tags.begin(), tags.end(), fit->info() ) != tags.end())
           cdt.delete_face(fit);
    } 
}

/** 
 * Based on CGAL output_to_medit, but for 2D meshes.
 * @see [output_to_medit] (https://github.com/CGAL/cgal/blob/master/Mesh_3/include/CGAL/IO/File_medit.h) 
 * @param os ostream of the output file  
 */
inline void Slice::output_slice_to_medit_(std::ostream& os)
{
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
  for( Vertex_iterator vit = cdt.vertices_begin(); vit != cdt.vertices_end();++vit)
  {
    V[vit] = inum++;
    vit->info() = 0; 
    os << *vit <<" "<< 0 <<'\n';
  }
  std::map<std::pair<Vertex_handle,Vertex_handle>,int> set_edges; 
  for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
  {
     for( int i =0; i<3; ++i)
     {  
         Edge eit(fit,i);
         Vertex_handle vh1 = fit->vertex(cdt.ccw(i));
         Vertex_handle vh2 = fit->vertex(cdt.cw(i));

         if (V[vh1] > V[vh2])
            set_edges[std::pair<Vertex_handle,Vertex_handle>(vh1,vh2)]=this->edges[eit];
         else
            set_edges[std::pair<Vertex_handle,Vertex_handle>(vh2,vh1)]=this->edges[eit];  
     }
  }
  os << "Edges\n" 
  << set_edges.size() << '\n';
  for(auto eit : set_edges) 
  {     
     os << V[eit.first.first] << " " << V[eit.first.second]  <<" " << eit.second <<std::endl;         
  }
  os << "Triangles\n" 
  << cdt.number_of_faces()<< '\n';
  for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
  {
        os << V[fit->vertex(0)]<< " " <<V[fit->vertex(1)]<<" " <<V[fit->vertex(2)]<< " " <<fit->info() <<"\n";
  }
  os << "End\n";
}

/** 
 * Constructor: 
 * Stores plane_3 and polylines to member variables. 
 * @param plane_3 3D plane that represents the slice    
 * @param polylines the constraints that determines the 2D mesh.
 */
inline Slice::Slice(Plane_3 plane_3,Polylines_2 &polylines) : plane(plane_3)
{
    min_sphere.add_polylines(polylines);
    constraints.insert( constraints.end(),polylines.begin(),polylines.end());
}

/** 
 * Adds constraints from another Slice object .
 * @param slice another Slice object
 */
inline void Slice::add_constraints(Slice &slice) 
{  
  add_constraints(slice.get_constraints());
}

/** 
 * Adds polyline to constraints .
 * @param polyline vector of points in 2D.
 */
inline void Slice::add_constraint(Polyline_2 &polyline) 
{  
     min_sphere.add_polyline(polyline);
     constraints.push_back( polyline );
}

/** 
 * Adds polylines to constraints .
 * @param polylines vector of vector of points in 2D.
 */
inline void Slice::add_constraints(Polylines_2 &polylines) 
{  
     min_sphere.add_polylines(polylines);
     constraints.insert( constraints.end(),polylines.begin(),polylines.end());
}

/** 
 * Adds constraints to the CGAL triangulation object cdt,
 * @param none
 */
inline void Slice::set_constraints() 
{    
     for (auto pol: constraints){cdt.insert_constraint(pol.begin(), pol.end());}
}

/** 
 * Creates 2D mesh.
 * @param mesh_resolution a double value divisor used with the  
          minimum bounding radius to determine the maximum  edge size.  
 */
inline void Slice::create_mesh(double mesh_resolution) 
{
     set_constraints(); 
     double r = min_sphere.get_bounding_sphere_radius();
     double longest_edge = r/mesh_resolution;
     Mesher mesher(cdt);    
     mesher.set_criteria(Criteria(0.125, longest_edge));

     std::cout << "Start  meshing."  << std::endl;
     mesher.refine_mesh(); 	
     std::cout << "Done  meshing."   << std::endl;
    
     for(CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
     {
         fit->info()=0;
         for( int i =0; i<3; ++i)
         {
            fit->vertex(i)->info()=0;
         }
     }
}

/** 
 * Slice a number of surfaces with the same plane, and store the constraints.
 * @param surfaces a vector of SVMTK Surface objects defined in Surface.h 
 */
template<typename Surface > 
void Slice::slice_surfaces(std::vector<Surface> surfaces) 
{
   for ( auto surf : surfaces ) 
   {
         std::shared_ptr<Slice> temp = surf.template mesh_slice<Slice>(this->plane);     
         this->add_constraints(*temp.get()); 
   }
}  

/** 
 * Transforms the 2D mesh to 3D surface mesh stores in a SVMTK Surface object.
 * @return surf SVMTK Surface object defined in Surface.h 
 */
template<typename Surface> 
std::shared_ptr<Surface> Slice::as_surface() 
{
  typedef std::vector<std::size_t> Face;
  typedef CDT::Vertex_handle Vertex_handle;

  std::vector<Point_3> points;
  std::vector<Face> faces;
  std::map<Vertex_handle, int> index_of_vertex;
  int i = 0;
  for(CDT::Point_iterator it = cdt.points_begin();  it != cdt.points_end(); ++it, ++i)
  {
       points.push_back(plane.to_3d(*it));  
       index_of_vertex[it.base()] = i;
  }
  for(CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
  {
       Face temp;
       temp.push_back( index_of_vertex[fit->vertex(0)] );
       temp.push_back( index_of_vertex[fit->vertex(1)] );
       temp.push_back( index_of_vertex[fit->vertex(2)] );
       faces.push_back(temp);
  }
   std::shared_ptr<Surface> surf(new  Surface(points, faces)); 
   return surf;  
}

/** 
 * Add tags to the facets in the 2D mesh based on overlapping surfaces and DefaultMap.                                      
 * @param surfaces a vector of SVMTK surface objects 
 * @return void
 * @see SubdomainMap.h 
 */
template<typename Surface> 
void Slice::add_surface_domains(std::vector<Surface> surfaces)
{
   if ( this->cdt.number_of_faces()==0)
   {
      std::cout<<"create mesh first"<< std::endl; 
      return;
   }
   DefaultMap map =DefaultMap();
   add_surface_domains(surfaces,map);
}

/** 
 * Add tags to the facets in the 2D mesh based on overlapping surfaces and SubdomainMaps.                                      
 * @param surfaces a vector of SVMTK surface objects 
 * @param map derived from SVMTK AbstractMap objects, @see SubdomainMap.h 
 * @return 
 */
template<typename Surface> 
void Slice::add_surface_domains(std::vector<Surface> surfaces, AbstractMap& map) 
{
   if ( this->cdt.number_of_faces()==0)
   {
      std::cout<<"create mesh first"<< std::endl; 
      return;
   }
   typedef boost::dynamic_bitset<>   Bmask;
   typedef std::pair<int,int>  Pid; 
   typedef std::map<Pid,int> Pid_map;
   int fn,fi;
   Pid_map pid_map_;
   Pid spp;
   std::map<Face_handle,Bmask> masks;
   int index_counter=1;
   
   for ( auto surf :  surfaces) 
   {
       typename Surface::Inside inside(surf.get_mesh());
       for(CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
       {
          Point_2 p2 =  CGAL::centroid(fit->vertex(0)->point() ,fit->vertex(1)->point(),fit->vertex(2)->point());  
          Point_3 p3 =  this->plane.to_3d(p2);
          CGAL::Bounded_side res = inside(p3);
          if (res == CGAL::ON_BOUNDED_SIDE)
             masks[fit].push_back(1);
          else
             masks[fit].push_back(0);
       }
   } 
   for (auto bit : masks)
   {
       bit.first->info()=map.index(bit.second);
       if (bit.first->info()==0)
          cdt.delete_face(bit.first);
   }
   Edge ei;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
      fi = fit->info(); 
      for( int i =0;i<3;++i)
      {
         Edge ei(fit,i); 
         if (fit->neighbor(i)->is_in_domain())
            fn= fit->neighbor(i)->info();
         else 
            fn=0;  
           
         if (fn!=fi)
         {     
            if (fn>fi)
                spp={fn,fi};
            else
                spp={fi,fn};
            std ::pair<Pid_map::iterator, bool> is_insert_successful = pid_map_.insert(std::make_pair(spp,index_counter));

            if(is_insert_successful.second)
            { 
                index_counter++;
            }
            this->edges[ei] = pid_map_[spp];    
        }
        else
        {
           this->edges[ei] =  0;
        }          
     }
  }        
}

/** 
 * Simplify constraints to a specific point density.
 * @param point_density the number of points pr. length
 * @return void. 
 */
inline void Slice::simplify(const double point_density )
{       
    Polylines_2 result;
    for ( auto c = this->constraints.begin(); c !=this->constraints.end(); ++c ) 
    {
          Polyline_2 temp; 
          double length = length_polyline(c->begin(),c->end()); 
          double adjustment = point_density*length/(double)(c->size());    
          CGAL::Polyline_simplification_2::simplify(c->begin(), c->end(), Cost() , Stop(adjustment), std::back_inserter(temp));
          result.push_back( temp ); 
    }
    this->constraints.clear();
    this->constraints=result;
}

/** 
 * Calculates and keeps the largest connected component.                                          
 * @param none 
 * @return void removes other connected components from stored mesh.  
 */
inline void Slice::keep_largest_connected_component() 
{
   std::vector<Face_handle> queue;
   std::map<Face_handle,bool> handled; 
   std::vector<std::vector<Face_handle>> connected_components;
   int num_cc=0;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
       handled[fit] = false; 
   }
   Face_handle fiq,fin;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
      std::vector<Face_handle> connnected_component;
      if (handled[fit]) continue;
      queue.push_back(fit);       
 
      while( !queue.empty()  )
      {
         fiq = queue.back(); 
         queue.pop_back();         
         if (handled[fiq]) continue;  
         handled[fiq] = true;
         connnected_component.push_back(fiq); 
         for( int i =0; i < 3; i++) 
         {        
            fin = fiq->neighbor(i); 
            if (  handled.find(fin)==handled.end() ) continue ;     
            queue.push_back(fin);
         }
      }
      num_cc++;
      connected_components.push_back(connnected_component); 
    }
    if (num_cc<2) return;
    std::sort(connected_components.begin(), connected_components.end(), sort_vectors_by_size());
    for ( auto ccit  = connected_components.begin()+1; ccit !=connected_components.end(); ccit++)
    {
          for ( auto fit = ccit->begin(); fit!=ccit->end(); fit++)
          {
              cdt.delete_face(*fit);

          }
    } 
}

/** 
 * Calculates and returns the number of connected components.                                         
 * @param none 
 * @return num_cc number of connected components. 
 */
inline int Slice::connected_components() 
{
   std::vector<Face_handle> queue;
   std::map<Face_handle,bool> handled; 
   int num_cc=0;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
       handled[fit] = false; 
   }
   Face_handle fiq,fin;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
      if (handled[fit]) continue;
      queue.push_back(fit);     
      while( !queue.empty()  )
      {
         fiq = queue.back(); 
         queue.pop_back();         
         if (handled[fiq]) continue;  
         handled[fiq] = true;
         for( int i =0; i < 3; i++) 
         {        
            fin = fiq->neighbor(i); 
            if (  handled.find(fin)==handled.end() ) continue;     
            queue.push_back(fin);
         }
      }
     num_cc++;      
   }
   return num_cc;
}

/** 
 * Saves the 2D mesh to file. 
 * Valid format are: off, stl, vtu and mesh(with tags)                                          
 * @param filename where the 2D mesh is to be stored
 */
inline void Slice::save(std::string outpath)
{
     std::string extension = outpath.substr(outpath.find_last_of(".")+1);
     if ( cdt.number_of_faces()==0 ) 
     {
        std::cout <<"The resulting mesh has no facet, and will not be saved"<< std::endl;
        return;         
     }

     if ( extension=="off")
     {
       std::ofstream out(outpath);
       CGAL::export_triangulation_2_to_off(out,cdt);
     } 
     else if ( extension=="stl")
     {
        write_STL(outpath);
     }
     else if ( extension=="vtu")
     {
        std::ofstream out(outpath);
        CGAL::write_vtu(out,cdt);  
     }
     else if ( extension=="mesh")
     {
         std::ofstream out(outpath);
         output_slice_to_medit_(out);
     }
}

/** 
 * Saves the 2D mesh to file in the stl format 
 * @param filename where the 2D mesh is to be stored
 */
inline void Slice::write_STL(const std::string filename)
{
    std::ofstream file(filename);
    file.precision(6);

    file << "solid "<< filename << std::endl;
    for (CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit )
    {
       Point_3 p1 = plane.to_3d( fit->vertex(0)->point() );
       Point_3 p2 = plane.to_3d( fit->vertex(1)->point() );
       Point_3 p3 = plane.to_3d( fit->vertex(2)->point() );
       
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






#endif 
