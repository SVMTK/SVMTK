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

//#include <CGAL/Alpha_shape_2.h>
//#include <CGAL/Alpha_shape_vertex_base_2.h>
//#include <CGAL/Alpha_shape_face_base_2.h>

#include <CGAL/utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/assertions.h>

#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>
// TODO: REWRITE
#include <algorithm> 
#include <iterator>
#include <queue>
#include <assert.h>  
#include <iterator>
//#include <vector>
//#include <set>
#include <fstream>

class Surface;

class Domain;


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




// TODO: move to new file
template< typename Kernel>
struct Minimum_sphere_2
{
   
    typedef typename CGAL::Min_sphere_of_spheres_d_traits_2< Kernel,typename Kernel::FT> Traits;
    typedef typename CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
    typedef typename Traits::Sphere                    Sphere;
 
    typedef typename std::vector<typename Kernel::Point_2> Polyline_2;
    typedef typename std::vector<Polyline_2> Polylines_2;  
  

    void add_polyline( const Polyline_2 &polyline) 
    {
         for (auto it=polyline.begin();it != polyline.end(); ++it)
         {
             S.push_back(Sphere(*it, 0.0));
         }
    }

    void add_polylines(const Polylines_2 &polylines)
    {
         for (auto it=polylines.begin();it != polylines.end(); ++it)
         {
             for ( auto pit=it->begin() ; pit!=it->end(); ++pit)
             { 
                 S.push_back(Sphere(*pit, 0.0));
             }
         }
     } 

     double get_bounding_sphere_radius()
     {
          Min_sphere ms(S.begin(), S.end());
           return CGAL::to_double(ms.radius());
     }
     private:
         std::vector<Sphere> S;
};



class Slice
{

    public :
     
       typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
       typedef Kernel::Point_2 Point_2;
       typedef Kernel::Line_2 Line_2;
       typedef Kernel::Point_3 Point_3;
       typedef Kernel::Segment_2 Segment;

       typedef CGAL::Triangulation_vertex_base_with_info_2<int,Kernel> Vb;
       typedef CGAL::Triangulation_face_base_with_info_2<int,Kernel> Fb_w_i;
       typedef CGAL::Constrained_triangulation_face_base_2<Kernel,Fb_w_i> C_fb_w_i;

       typedef CGAL::Delaunay_mesh_face_base_2<Kernel,C_fb_w_i> Fb;

       typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
       typedef Tds::Vertex_circulator Vertex_circulator;
       typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds,CGAL::Exact_predicates_tag> CDT1; //??
       typedef CGAL::Constrained_triangulation_plus_2<CDT1> CDT;
       
       typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
       typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;
       typedef Kernel::FT FT;
       typedef std::vector<Point_2> Polyline_2;
       typedef std::vector<Polyline_2> Polylines_2;      

       typedef Kernel::Plane_3 Plane_3; 
       typedef CGAL::Polygon_2<Kernel> Polygon_2;
       typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_wh2;
       typedef std::vector<Polygon_2> Polygons_2;
       //typedef CGAL::Polyline_simplification_2::Stop_above_cost_threshold Stop;
       typedef CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold Stop;
       typedef CGAL::Polyline_simplification_2::Squared_distance_cost Cost;
       typedef CDT::Face_handle Face_handle ;
       typedef CDT::Vertex_handle Vertex_handle ;
       typedef CDT::Vertex_iterator Vertex_iterator;
       typedef CDT::Face_iterator Face_iterator;
       typedef CDT::Edge Edge;
       // TODO: prequisit of slice = Plane_3
       Slice(Plane_3 plane ,Polylines_2 &polylines) ;
       Slice(Plane_3 plane_3) : plane(plane_3) {};
       Slice(){}
       ~Slice(){} 


       void set_plane(Plane_3 inplane){ this->plane = inplane;}
       Plane_3& get_plane(){return this->plane;}

       void add_surface_domains(std::vector<Surface> surfaces, AbstractMap& map) ; //todo rename
       void add_surface_domains(std::vector<Surface> surfaces) ; //todo rename
       void slice_surfaces(std::vector<Surface> surfaces) ;

       void remove_subdomain(std::vector<int> tags); 
       void create_mesh(double mesh_resolution);               
       void simplify(double stop_crit); 
       
    

       int connected_components(); 
       void keep_largest_connected_component(); 

       void add_constraints(Polylines_2 &polylines, bool hole=false); 
       void add_constraints(Slice &slice, bool hole=false);
       void add_constraint(Polyline_2 &polyline, bool hole=false);



       void repair_polylines(Polylines_2& polylines_bad) ;

       void clear_costraints(){ constraints.clear(); return; } 
       int  num_constraints() { return constraints.size();}
       void set_constraints() ;
       void repair_costraints() { Polylines_2 temp = constraints; clear_costraints() ;std::sort( temp.begin(), temp.end(), sort_vectors_by_size() ); repair_polylines(temp);}


       std::map<Edge,int>& get_edges(){this->edges;}  
       Polylines_2& get_constraints() { return constraints;}

       std::shared_ptr<Surface> as_surface(); //Forward called
       void save(std::string outpath);
       void output_slice_to_medit_(std::ostream& os);
       void write_STL(const std::string filename);
       void simplify_polylines(Polylines_2& polylines, double point_density=0.4 );
       // Structs 
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
       std::map<Edge,int> edges; // after failing to use info.

};

void Slice::remove_subdomain(std::vector<int> tags) 
{
    for(CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
    {
        if(std::find(tags.begin(), tags.end(), fit->info() ) != tags.end())
           cdt.delete_face(fit);
    } 
}

void Slice::output_slice_to_medit_(std::ostream& os)
{
  // Based on CGAL output_to_medit, but for 2D meshes
 //TODO ; improve

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

  int fi, fn;
  int edge_counter;
  for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
  {
     for( int i =0 ; i<3 ; ++i)
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


Slice::Slice(Plane_3 plane_3,Polylines_2 &polylines) : plane(plane_3)
{
    min_sphere.add_polylines(polylines);
    simplify_polylines(polylines);
    constraints.insert( constraints.end(),polylines.begin(),polylines.end());
}



void Slice::add_constraints(Slice &slice, bool hole) 
{  
  add_constraints(slice.get_constraints(), hole);
}


void Slice::add_constraint(Polyline_2 &polyline, bool hole) 
{  
     min_sphere.add_polyline(polyline);
     if (hole)
     {
        if  ( static_cast<int>(CGAL::orientation_2(polyline.begin(),polyline.end()))== -1) 
        {
           std::reverse(polyline.begin(),polyline.end());
        }
        Point_2 c2 = CGAL::centroid(polyline.begin(), polyline.end(), CGAL::Dimension_tag<0>());
        
        constraints.push_back(polyline);
                     
     }
     else 
     {
         constraints.push_back( polyline );
     }
}
void Slice::add_constraints(Polylines_2 &polylines, bool hole) 
{  
     min_sphere.add_polylines(polylines);
     constraints.insert( constraints.end(),polylines.begin(),polylines.end());
}

void Slice::simplify( double stop_crit )
{
     if(constraints.size() < 1 or constraints[0].size() < 10)
     {
       return;
     }

     Polylines_2 temp;
     for ( auto  pol = constraints.begin(); pol != constraints.end(); ++pol ) 
     {
         Polyline_2 result;
         CGAL::Polyline_simplification_2::simplify(pol->begin(), pol->end(), Cost() , Stop(stop_crit), std::back_inserter(result));
         temp.push_back(result);
  
     }
     constraints = temp;
}

void Slice::set_constraints() 
{    for (auto pol: constraints){cdt.insert_constraint(pol.begin(), pol.end());}
}

void Slice::create_mesh(double mesh_resolution) 
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
         for( int i =0 ; i<3 ; ++i)
         {
            fit->vertex(i)->info()=0;
         }
     }


}


void Slice::simplify_polylines(Polylines_2& polylines, double point_density )
{
   // ----------------------------------------
   // In general, we want the that the density points/length to be less than 0.5, which gives 
   // the edge size of 10 point to be greater than 50 length units
   // -----------------------------------------             
    Polylines_2 result;
    for ( auto c = polylines.begin(); c !=polylines.end() ; ++c ) 
    {
          Polyline_2 temp; 
          double length = length_polyline(c->begin(),c->end()) ; 
          double adjustment = point_density*length/(double)(c->size());    
          CGAL::Polyline_simplification_2::simplify(c->begin(), c->end(), Cost() , Stop(adjustment), std::back_inserter(temp));
          result.push_back( temp ) ; 
    }
    polylines.clear();
    polylines=result;


}

void Slice::keep_largest_connected_component() 
{
   std::vector<Face_handle> queue;
   std::map<Face_handle,bool> handled; 
   std::vector<std::vector<Face_handle>> connected_components;
   int num_cc=0;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
       handled[fit] = false ; 
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
         if (handled[fiq]) continue ;  
         handled[fiq] = true;
         connnected_component.push_back(fiq); 
         for( int i =0 ; i < 3 ; i++) 
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
          for ( auto fit = ccit->begin() ; fit!=ccit->end(); fit++)
          {
              cdt.delete_face(*fit);

          }

    } 

}



int Slice::connected_components() 
{
   std::vector<Face_handle> queue;
   std::map<Face_handle,bool> handled; 
   int num_cc=0;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
       handled[fit] = false ; 
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
         if (handled[fiq]) continue ;  
         handled[fiq] = true;
         for( int i =0 ; i < 3 ; i++) 
         {        
            fin = fiq->neighbor(i); 
            if (  handled.find(fin)==handled.end() ) continue ;     
            queue.push_back(fin);
         }
      }
     num_cc++;
      
   }
   return num_cc;
}




void Slice::save(std::string outpath)
{
     std::string extension = outpath.substr(outpath.find_last_of(".")+1);
   
     if ( cdt.number_of_faces()==0 ) 
     {
        std::cout <<"The resulting mesh has no facet, and will not be saved"<< std::endl;
        return ;         
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
         output_slice_to_medit_(out) ;
     }


}


void Slice::write_STL(const std::string filename)
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
