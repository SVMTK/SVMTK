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
#include <boost/dynamic_bitset.hpp>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/IO/Triangulation_off_ostream_2.h>
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
#include <assert.h>  
#include <iterator>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_vtu.h>
#include <algorithm> 

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Polygon_with_holes_2.h>

// TODO: REWRITE
#include <iterator>
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
       Slice(Plane_3 plane_3) : plane(plane_3) {};
       Slice(){}

       ~Slice(){}
 
       Slice(Plane_3 plane ,Polylines_2 &polylines) ;
       void add_surface_domains(std::vector<Surface> surfaces, AbstractMap& map) ; //todo rename
       void add_surface_domains(std::vector<Surface> surfaces) ; //todo rename
       void slice_surfaces(std::vector<Surface> surfaces) ;
       void remove_subdomain(std::vector<int> tags); 
       void write_STL(const std::string filename);
       void create_mesh(double mesh_resolution);               
       void simplify(double stop_crit); 
       void set_plane(Plane_3 inplane){ this->plane = inplane;}
       Plane_3& getPlane(){return this->plane;}
       void add_subdomains(Domain& domain);
       void save(std::string outpath);

       int connected_components(); 
       void keep_largest_connected_components(); 
       void remove_isolated_vertices();

       void add_constraints(Polylines_2 &polylines, bool hole=false); 
       void add_constraints(Slice &slice, bool hole=false);
       void add_constraint(Polyline_2 &polyline, bool hole=false);
       void repair_polylines(Polylines_2& polylines_bad) ;
       void set_constraints() ;
       std::map<Edge,int>& get_edges(){this->edges;}
       void repair_costraints() { Polylines_2 temp = constraints; clear_costraints() ;std::sort( temp.begin(), temp.end(),compare_size() ); repair_polylines(temp);}
       void clear_costraints(){constraints.clear();} 
     
       Polylines_2& get_constraints() { return constraints;}

       std::shared_ptr<Surface> export_3D(); //Forward called
       void output_slice_to_medit_(std::ostream& os);
       int  num_constraints() { return constraints.size();}
       void simplify_polylines(Polylines_2& polylines, double point_density=0.4 );
       struct polyline_endpoints{

                  polyline_endpoints( const std::vector<Point_2> & in ) : current(in) {}
                  bool operator()(const std::vector<Point_2> & a, const std::vector<Point_2>& b)
                       { 
                          double tmp1 = static_cast<double>(CGAL::min(CGAL::squared_distance(a.front(),current.back()),CGAL::squared_distance(a.back(),current.back() ) ));
                          double tmp2 = static_cast<double>(CGAL::min(CGAL::squared_distance(b.front(),current.back()),CGAL::squared_distance(b.back(),current.back() ) ));

                          return tmp1 < tmp2; 
                       }
                  private:
                      std::vector<Point_2> current;
                   
       };

       struct search_knot{
                  search_knot(const std::vector<Point_2> & b) : current(b) {} 
                  bool operator()(const std::vector<Point_2> & a)
                       { 
                          if ( a.front()== current.front() and  a.back() == current.back()  ) 
                              return true;

                          return false; 
                       }
                  private:
                      std::vector<Point_2> current;
                   
       };
       struct search_lens{
                       search_lens(const std::vector<Point_2> & b) : current(b) {} 
                       bool operator()(const std::vector<Point_2> & a)
                       { 
                          if ( a.front()== current.front() and  a.back() == current.back()  ) 
                              return true;
                          else if ( a.front()== current.back() and  a.back() == current.front()  ) 
                              return true;
                          return false; 
                       }
                  private:
                      std::vector<Point_2> current;
                   
       };
       struct compare_size {
                  bool operator()(const std::vector<Point_2> & a, const std::vector<Point_2> & b)
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
  for( Vertex_iterator vit = cdt.vertices_begin();
       vit != cdt.vertices_end();
       ++vit)
  {
    V[vit] = inum++;
    vit->info() = 0; 
    os << *vit <<" "<< 0 <<'\n';
  }
  
  std::map<std::pair<Vertex_handle,Vertex_handle>,int> set_edges; 

  int fi, fn;
  int id_edge;
  int edge_counter;
  // TODO : rename ib to fit 
  id_edge=0;
  for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
  {
     for( int i =0 ; i<3 ; ++i)
     {  
         Edge  eit(fit,i);
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
  for( auto eit : set_edges) 
  {     
         os << V[eit.first.first] << " " << V[eit.first.second]  <<" " << eit.second <<std::endl;         
  }
  os << "Triangles\n" 
  << cdt.number_of_faces()<< '\n';
  for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
  {
     
        os << V[fit->vertex(0)]<< " " <<V[fit->vertex(1)]<<" " <<V[fit->vertex(2)]<< " " <<fit->info() <<"\n";
   
  }





  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End\n";

}
// TODO : slice(Plane_3 () ) , std::vector<Surface> surfaces, AbstractMap& map) 
// slice.slice_surfaces();
// slice.slice_mesh() ;
/*

void Slice::slice_mesh(Domain domain)
{

   std::vector<Domain::Point_3> points;
   std::vector<Domain::Face> faces;

   for ( auto i : domain.number_of_patches() ) // 
   {
      points.clear(); 
      faces.clear();
      facets_in_complex_3_to_triangle_soup_(domain.get_mesh(),Domain::Surface_patch_index(i.first,i.second),points,faces);
      std::shared_ptr<Surface> surf( new Surface(points,faces)); 
      std::shared_ptr<Slice> temp = surf->mesh_slice(this->plane);     // nerror here   
      this->add_constraints(*temp.get()); 
    
   }
   this->add_subdomains(domain);
}  
  



void Slice::add_subdomains(std::vector<Surface> surfaces)
{
     DefaultMap* map = new DefaultMap();
     add_subdomains(surfaces, map); 

} */


// TODO : slice( Domain domain)   

Slice::Slice(Plane_3 plane_3,Polylines_2 &polylines) : plane(plane_3)
{
   
    min_sphere.add_polylines(polylines);
    simplify_polylines(polylines);
    //std::sort(polylines.begin(), polylines.end(),compare_size() );

    //repair_polylines(polylines);
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
{   

    
//     repair_costraints();
     
     
     for ( auto  pol : constraints ) 
     {         
               cdt.insert_constraint(pol.begin(), pol.end()); 
               
     }
}



void Slice::create_mesh(double mesh_resolution) 
{
     std::cout << "Setting constraints  " << constraints.size() << std::endl;
     set_constraints(); 
     std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
     double r = min_sphere.get_bounding_sphere_radius();
     double longest_edge = r/mesh_resolution;
     
     Mesher mesher(cdt);    
     

     mesher.set_criteria(Criteria(0.125, longest_edge));

     std::cout << "Start  meshing "  << std::endl;
     mesher.refine_mesh(); 	

     std::cout << "Done  meshing   " << cdt.number_of_faces() << std::endl;


      
     for(CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
     {
         fit->info()=0;
         for( int i =0 ; i<3 ; ++i)
         {
            fit->vertex(i)->info()=0;
         }
 


     }


}
template< typename Polyline > // TODO:change to iterator 
double length_polyline( Polyline& polyline)
{
  double length = 0.0;
  for (auto cit = polyline.begin() ; cit!=polyline.end(); ++cit)
  {
     length += static_cast<double>( CGAL::sqrt(CGAL::squared_distance(*cit, *(cit+1) ) ) ); 

  } 
        
  return length; 


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
          double length = length_polyline(*c) ; 
          double adjustment = point_density*length/(double)(c->size());  
      
          CGAL::Polyline_simplification_2::simplify(c->begin(), c->end(), Cost() , Stop(adjustment), std::back_inserter(temp));
          result.push_back( temp ) ; 
    }
    polylines.clear();
    polylines=result;


}


void Slice::repair_polylines(Polylines_2& polylines_bad) 
{


  
    // TODO : rename cyclic and acyclic to open and closed ??
    Polylines_2 cyclic;
    Polylines_2 acyclic;
    Polylines_2 lp;
    int step=0;
    Polyline_2 pline;
    std::set<Polylines_2::iterator> toremove;
    std::cout << "step"<< step++  << std::endl;
    for ( auto c = polylines_bad.begin(); c !=polylines_bad.end() ; ++c ) 
    {

       if (c->size() <  3) 
           continue;

       if ( c->front() == c->back()) 
       {  
           cyclic.push_back(*c) ;
       }
       else  
       {
          acyclic.push_back(*c);
       }         
    }
    std::cout << "step"<< step++  << std::endl;
    // ----------------------------------------
    //  Find pairs of open Polylines that create a lense object 
    //  and add to closed polygons instead.
    // -----------------------------------------
    for ( int i = 1 ; i < acyclic.size() ; ++i ) 
    {
       auto it  = std::find_if(acyclic.begin()+i, acyclic.end(), search_lens(acyclic[i-1]) ); 
       if ( it!= acyclic.end() ) 
       { 
            Polyline_2 temp;
            toremove.insert(it); 
            toremove.insert(acyclic.begin()+i-1);   
  
            temp.insert(temp.end(), acyclic[i-1].begin(),acyclic[i-1].end()-1);
            temp.insert(temp.end(), it->begin(), it->end());
            cyclic.push_back(temp);
       }    
    }
    std::cout << "step"<< step++  << std::endl;    
    for ( auto c = toremove.begin(); c!= toremove.end(); ++c) 
    {
         acyclic.erase(*c);
    }
    // ----------------------------------------
    // Finds open polygons matching endpoints to  form a closed polygon of the boundary. TODO: matching or least distance away
    // Polylines that causes a 1D constraint on the boundary are removed 
    // Polylines that are not in the closed polygon are removed.
    // -----------------------------------------
    int ac_size = acyclic.size();
    std::cout << "step"<< step++  << acyclic.size() <<std::endl;


     
    for ( int i = 1 ; i < ac_size ; ++i ) 
    {
        std::sort(acyclic.begin()+i, acyclic.end(),polyline_endpoints(acyclic[i-1]) );
        
        if ( CGAL::squared_distance(acyclic[i-1].back(),acyclic[i].back()) < CGAL::squared_distance(acyclic[i-1].back(),acyclic[i].front() )  )
        {
             std::reverse(acyclic[i].begin(),acyclic[i].end());
        }
        if ( acyclic[i-1].front() == acyclic[i].back() or acyclic[i-1].front() == acyclic[i].front())
        {

            acyclic.erase(acyclic.begin()+i-1) ;
            ac_size--;
            i--;



        }
        if ( acyclic[i].back()==acyclic[0].front() and i>1) 
        { 
            lp.insert(lp.end(),acyclic.begin()+i+1, acyclic.end());
            acyclic.erase(acyclic.begin()+i+1, acyclic.end()) ;
            break;
        }
    }
    std::cout << "step"<< step++  << std::endl;


    if ( acyclic.size()<3)
    {
        add_constraints(cyclic);
        return;
    }

    for ( auto c = acyclic.begin(); c !=acyclic.end() ; ++c ) 
    {

          Polyline_2 temp; 
          double length = length_polyline(*c) ; 
          double adjustment = 0.4*length/(double)(c->size());
         // ----------------------------------------
         // In general, we want the that the density points/length to be less than 0.5, which gives 
         // the edge size of 10 point to be greater than 50 length units
         // -----------------------------------------             
          CGAL::Polyline_simplification_2::simplify(c->begin(), c->end(), Cost() , Stop(adjustment), std::back_inserter(temp));
          pline.insert( pline.end(), c->begin(),c->end()-1 ) ; 
    }

    std::cout << "step"<< step++  << std::endl;

    if  ( static_cast<int>(CGAL::orientation_2(pline.begin(),pline.end()))== 1) 
          std::reverse(pline.begin(),pline.end());
 
    std::cout << "step"<< step++  << std::endl;
    // ----------------------------------------
    // Finds 
    // 
    // -----------------------------------------
    Polyline_2 result;
    CGAL::Polyline_simplification_2::simplify(pline.begin(), pline.end(), Cost() , Stop(0.8), std::back_inserter(result));


    constraints.insert(constraints.begin(),acyclic.begin(),acyclic.begin() ) ;
    
    add_constraints(cyclic);
    

    if ( lp.size() > 2 )  
    { 
         std::cout << "step"<< step++  << lp.size() << std::endl;
         std::sort(lp.begin(), lp.end(),compare_size());
         repair_polylines(lp);
    }


}
void Slice::remove_isolated_vertices() 
{

  // TODO FIX OR REMOVE
  std::map<Vertex_handle,bool> handled;
  typedef CDT::All_vertices_iterator All_vertices_iterator;
  typedef CDT::Finite_vertices_iterator Finite_vertices_iterator;
  for (Finite_vertices_iterator  vit = cdt.finite_vertices_begin() ; vit != cdt.finite_vertices_end() ; vit++)
  {
      handled[vit]=false;
  }
  for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
  {
     for( int i =0 ; i<3 ; ++i)
     {  
         Vertex_handle vh1 = fit->vertex(cdt.ccw(i));
         handled[vh1]=true;   
     }  
  }
  std::cout << cdt.number_of_vertices()<<std::endl;

  for (auto vit : handled) 
  {

     if (!vit.second && !cdt.is_infinite(vit.first) && vit.first!=Vertex_handle())  
     {
       std::cout << "HELP" << cdt.are_there_incident_constraints(vit.first) << cdt.is_infinite(vit.first) << std::endl;
       cdt.remove_incident_constraints(vit.first);
       std::cout << "HELP" << cdt.are_there_incident_constraints(vit.first) << cdt.is_infinite(vit.first) << std::endl;
       cdt.remove(vit.first) ;

     }
  }
  std::cout << cdt.number_of_vertices() << std::endl;

}

void Slice::keep_largest_connected_components() 
{
  int num_cc = connected_components();
  if (num_cc==1) return; 

   typedef CDT::Finite_vertices_iterator Finite_vertices_iterator;
   for (Finite_vertices_iterator  vit = cdt.finite_vertices_begin() ; vit != cdt.finite_vertices_end() ; vit++)
   {
     if (!vit->face()->is_in_domain()) continue; 
     
     if (vit->info()==0) cdt.delete_face(vit->face());
        
   }

}
int Slice::connected_components() 
{
   typedef CDT::Finite_vertices_iterator Finite_vertices_iterator;
   std::map<Face_handle,bool> handled;
   /*std::map<Vertex_handle,bool> handled;
   for (Finite_vertices_iterator  vit = cdt.finite_vertices_begin() ; vit != cdt.finite_vertices_end() ; vit++)
   {
      handled[vit]=false;
   }
      Vertex_circulator vc,vc0;

   std::vector<Vertex_handle> queue;
   */
   for(auto fit = cdt.all_faces_begin(); fit != cdt.all_faces_end(); ++fit) 
   {
       handled[fit]=false;
   }
   std::vector<Face_handle> queue;

   int id=1;
   int iter=0;
   int iterw=0;
   //typedef CDT::Finite_edges_iterator Finite_edges_iterator;
   typedef CDT::Vertex_circulator Vertex_circulator;
   typedef CDT::Finite_vertices_iterator Finite_vertices_iterator;
   Face_handle fic,fin;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
      std::cout << "IS";
      if (handled[fit]) continue;
      std::cout << "It";
      queue.push_back(fit);

      while( queue.size()>0  )
      {
         std::cout << "WORK";
         fic = queue.back();   
         std::cout << "ING";
         queue.pop_back();
         std::cout << "OR";
         if (handled[fic])
         { 
          
           continue;
         }
         handled[fic]=true;
         std::cout << "NOT";
         if(!fic->is_in_domain()) continue;
         std::cout << "?"<<std::endl;
         for( int i =0 ; i < 3 ; i++) 
         { 
            std::cout << "Y";
            fin = fic->neighbor(i);
            std::cout << "E";
         
            if (handled[fin]) continue;
            std::cout << "S";
            queue.push_back(fin);
            std::cout << "!" << std::endl;
         }
      std::cout <<"Queue" <<queue.size() << std::endl;
     }
     id++;
      
   }
  /* for (Finite_vertices_iterator  vit = cdt.finite_vertices_begin() ; vit != cdt.finite_vertices_end() ; vit++)
   {

    if (handled[vit]) continue;
    queue.push_back(vit);

    while( !queue.empty() )
    {
     Vertex_handle vic = queue.back();   
     queue.pop_back();

     if (handled[vic]) continue;

     handled[vic]=true;
     vic->info() = id;

     if(!vic->face()->is_in_domain()) continue ;

     vc0 = vc = cdt.incident_vertices(vic,vic->face());


     if ( !vc.is_empty() )
     { if(vc != NULL)
        { do
          { 
            if (!handled[vc])
                queue.push_back(vc);
             
            ++vc;           
          } while(vc!=vc0);
         }
       }
       std::cout <<"Queue" <<queue.size() << std::endl;
     }
     id++;
   }*/
   return --id;
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
