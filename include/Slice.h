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
#include <algorithm> 

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Polygon_with_holes_2.h>

class Surface;



// TODO: move to new file
template< typename Kernel>
struct Minimum_sphere_2
{
   
    typedef typename CGAL::Min_sphere_of_spheres_d_traits_2< Kernel,typename Kernel::FT> Traits;
    typedef typename CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
    typedef typename Traits::Sphere                    Sphere;
 
    typedef typename std::vector<typename Kernel::Point_2> Polyline_2;
    typedef typename std::vector<Polyline_2> Polylines_2;  
  
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
    // EXPERIMENTAL CLASS
    // TODO:CLEANUP
    public :
     
       typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
      
       typedef Kernel::Point_2 Point_2;
       typedef Kernel::Line_2 Line_2;
       typedef Kernel::Point_3 Point_3;
       typedef Kernel::Segment_2 Segment;
       typedef CGAL::Triangulation_vertex_base_2<Kernel> Vb;
       typedef CGAL::Delaunay_mesh_face_base_2<Kernel> Fb;
       typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
       typedef Tds::Vertex_circulator Vertex_circulator;
       typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, Tds,CGAL::Exact_predicates_tag> CDT; //??
       //typedef CGAL::Constrained_triangulation_plus_2<CDT1> CDT;
       
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

       typedef CGAL::Alpha_shape_2<CDT> Alpha_shape_2;
       Slice(){}
       ~Slice(){}
 
       Slice(Slice &slice){constraints = slice.get_constraints() ; } 

       Slice(const Polylines_2 &polylines) ;



       void write_STL(const std::string filename);
       void create_mesh(double mesh_resolution);               
       void simplify(double stop_crit); 
       void set_plane(Plane_3 inplane){ plane = inplane;}


       void save(std::string outpath);

       void add_constraints(Polylines_2 &polylines) {min_sphere.add_polylines(polylines); constraints.insert( constraints.end(),polylines.begin(),polylines.end());}

       void clear_costraints(){constraints.clear();} 

       Polylines_2& get_constraints() { return constraints;}
       Polygon_2& get_boundary() { return boundary;}

       void repair_domain(Polylines_2& polylines_bad , bool is_boundary); 
       std::shared_ptr<Surface> export_3D();

       int  num_constraints() { return constraints.size();}

       // TODO : make functions operational after usage of  polyline
       void add_constraints(Slice &slice, bool hole);
       void keep_component(int next );
       void add_holes(Polylines_2& closed_polylines );
       // TODO : remove redundant 
       //void find_holes(int min_num_edges);
       void remove_bad_constraints( int min_num_edges);
 


       // TODO: move structs to new file ? use template ?  
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
       Polyline_2 seeds;
       Polylines_2 constraints; // polygons_2D ?? 
       Polygon_2 boundary;  
       CDT cdt;
       Plane_3 plane;// Store the plane  equation -


};
//TODO :remove Polygon dependence


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



void Slice::repair_domain(Polylines_2& polylines_bad , bool is_boundary=false) 
{
    // TODO : rename cyclic and acyclic to open and closed ??
    Polylines_2 cyclic;
    Polylines_2 acyclic;
    Polylines_2 lp;

    Polyline_2 pline;
    std::set<Polylines_2::iterator> toremove;
    // ----------------------------------------
    //  Polylines are divived into  closed polygons and open polygons TODO:rename 
    //  Polylines less than 10 points are not considered TODO: make optional
    // -----------------------------------------
    for ( auto c = polylines_bad.begin(); c !=polylines_bad.end() ; ++c ) 
    {
       if (c->size() <  16) 
           break;

       if ( c->front() == c->back()) 
       {  
           cyclic.push_back(*c) ;
       }
       else  
       {
          acyclic.push_back(*c);
       }         
    }
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
    
    for ( auto c = toremove.begin(); c!= toremove.end(); ++c) 
    {
         acyclic.erase(*c);
    }
    // ----------------------------------------
    // Finds open polygons matching endpoints to  form a closed polygon of the boundary. TODO: mateching or least distance away
    // Polylines that causes a 1D constraint on the boundary are removed 
    // Polylines that are not in the closed polygon are removed.
    // -----------------------------------------
    int ac_size = acyclic.size();

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
    // ----------------------------------------
    // Combines open polygons to form a closed polygon of the boundary.
    // Density based 
    // -----------------------------------------    

    for ( auto c = cyclic.begin(); c !=cyclic.end() ; ++c ) 
    {
          Polyline_2 temp;
          double length = length_polyline(*c) ; 
          double adjustment = 0.4*length/(double)(c->size()); 
          if( (double)(c->size()) >length*0.5)  
          {
              CGAL::Polyline_simplification_2::simplify(c->begin(), c->end(), Cost() , Stop(adjustment), std::back_inserter(temp));
              *c=temp;
          }
          length = length_polyline(*c) ;
    }

    if ( acyclic.size()<3)
    {
        add_holes(cyclic);
        return;
    }
   
    //TODO ; remove pline and use acyclic  
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


    // ----------------------------------------
    //  Fixes the orientation and the number of points in the polygon
    // 
    // -----------------------------------------
    // TODO : remove



    if  ( static_cast<int>(CGAL::orientation_2(pline.begin(),pline.end()))== 1) 
          std::reverse(pline.begin(),pline.end());
 

    // ----------------------------------------
    // Finds 
    // 
    // -----------------------------------------
    Polyline_2 result;
    CGAL::Polyline_simplification_2::simplify(pline.begin(), pline.end(), Cost() , Stop(0.8), std::back_inserter(result));

    if (is_boundary) 
    { 
          //constraints.push_back(result)
          boundary.insert(boundary.vertices_end(),result.begin(), result.end()) ;
    }

    add_holes(cyclic);
    

    if ( lp.size() > 2 )  
    { 
         std::sort(lp.begin(), lp.end(),compare_size());
         repair_domain(lp,false);
    }

}

Slice::Slice(const Polylines_2 &polylines) 
{

    min_sphere.add_polylines(polylines);

    Polylines_2 temp; 
    temp = polylines;
   
    std::sort( temp.begin(), temp.end(),compare_size() );
     
    if ( temp[0].size() > 1 ) 
    { 

       if (temp[0].front() != temp[0].back() ) 
       {
           repair_domain(temp,true);
       }
       else 
       { 
  
           boundary.insert(boundary.vertices_end(),temp[0].begin(), temp[0].end() ) ;
           temp.erase(temp.begin());
           repair_domain(temp,false);
       }

    }

    
}


std::shared_ptr<Surface> Slice::export_3D() 
{
  // FIXME
  /*typedef std::vector<std::size_t> Face ;
  typedef CDT::Vertex_handle Vertex_handle;

  std::vector<Point_3> points ;
  std::vector<Face> faces;
  std::map<Vertex_handle, int> index_of_vertex;
  int i = 0;
  for(CDT::Point_iterator it = cdt.points_begin() ;  it != cdt.points_end(); ++it, ++i)
  {
       points.push_back(plane.to_3d(*it));  
       index_of_vertex[it.base()] = i;
  }
  
  for(CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
  {
       Face temp;
       temp.push_back( index_of_vertex[fit->vertex(0)] ) ;
       temp.push_back( index_of_vertex[fit->vertex(1)] ) ;
       temp.push_back( index_of_vertex[fit->vertex(2)] ) ;
       faces.push_back(temp);
  }
   std::shared_ptr<Surface> surf(new  Surface(points, faces)) ;
   return surf;  */
}

void Slice::keep_component(int next )
{ 
    // FIXME
    if ( next ==0)
    {
       constraints.clear();
    }
    else 
    {
        if ( constraints.size()-1 < next) 
        {
            next = constraints.size()-1; 
        }
     
        Polyline_2 temp= constraints[next-1]; 
        boundary.clear();
        boundary.insert(boundary.vertices_end(),temp.begin(), temp.end()) ;
    } 

}



void Slice::add_holes(Polylines_2& closed_polylines )
{   
    //FIXME  
    if ( boundary.is_empty() )
         return; 
     for ( auto pol =  closed_polylines.begin(); pol != closed_polylines.end(); ++pol)
     {
            Polyline_2 temp = *pol;
            if  ( static_cast<int>(CGAL::orientation_2(temp.begin(),temp.end()))== -1) 
            {
                  std::reverse(temp.begin(),temp.end());
            }
            Point_2 c2 = CGAL::centroid(temp.begin(), temp.end(), CGAL::Dimension_tag<0>());
        
            if ( CGAL::bounded_side_2(temp.begin(), temp.end(), c2 , Kernel())==CGAL::ON_BOUNDED_SIDE and  boundary.has_on_bounded_side(c2) )        
            {
                constraints.push_back(temp);
                seeds.push_back(c2); 
            }
     }
    
} 

void Slice::add_constraints(Slice &slice, bool hole=false) 
{  
     if (hole)
     {
           for ( auto pol =  slice.get_constraints().begin(); pol != slice.get_constraints().end(); ++pol)
           {
                Point_2 c2 = CGAL::centroid(pol->begin(), pol->end(), CGAL::Dimension_tag<0>());
                seeds.push_back(c2); 
           }
     }
     add_constraints(slice.get_constraints());
}


void Slice::remove_bad_constraints( int min_num_edges)
{
      if(constraints.size() < 2 or constraints[0].size() < 20)
      {
        constraints.clear();
        return;
      }

      for ( auto  pol = std::next(constraints.begin())  ; pol != constraints.end();) 
      {   
         Point_2 c2 = CGAL::centroid(pol->begin(), pol->end(), CGAL::Dimension_tag<0>());

         if( !CGAL::is_simple_2( pol->begin(), pol->end(), Kernel()) or pol->size() < min_num_edges)
         {
            pol = constraints.erase(pol);
         }
         else 
         {
           ++pol;
         } 

      }


} 

//TODO: rmeove
/*void Slice::find_holes( int min_num_edges)
{

     if(constraints.size() < 2)
     {
      //std::cout  << "There are no holes."<<  std::endl;
      return;
     }
     
     Polyline_2 temp = constraints[0];
     int j=0;
     Point_2 c2 = CGAL::centroid(temp.begin(), temp.end(), CGAL::Dimension_tag<0>());
     if ( CGAL::bounded_side_2(temp.begin(), temp.end(), c2 , Kernel())==CGAL::ON_UNBOUNDED_SIDE)
     {
       std::cout << "Bad slice" << std::endl;
     } 


     for ( auto  pol = std::next(constraints.begin())  ; pol != constraints.end();) 
     {
         Point_2 c2 = CGAL::centroid(pol->begin(), pol->end(), CGAL::Dimension_tag<0>());
         if ( CGAL::bounded_side_2(temp.begin(), temp.end(), c2 , Kernel())==CGAL::ON_BOUNDED_SIDE  and                                // inside the largest polyline 
              CGAL::bounded_side_2(pol->begin(), pol->end(), c2 , Kernel())==CGAL::ON_BOUNDED_SIDE  and  pol->size() > min_num_edges ) // inside its own polyline i.e. have enclosed area          
         {
            seeds.push_back(c2);
            ++pol;
         } 
         else if ( CGAL::bounded_side_2(temp.begin(), temp.end(), c2 , Kernel())==CGAL::ON_UNBOUNDED_SIDE or pol->size() < min_num_edges)           
         {
            pol = constraints.erase(pol);
        
         }
         else
         {
            ++pol;
         }
     }

   
}*/


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



void Slice::create_mesh(double mesh_resolution) 
{
    // 

    if ( boundary.is_empty() )
         return; 


    cdt.insert_constraint(boundary.vertices_begin() , boundary.vertices_end(), true ) ;

     for ( auto  pol = constraints.begin(); pol != constraints.end(); ++pol ) 
     {         
               cdt.insert_constraint(pol->begin(), pol->end(),true); //
     }

     double r = min_sphere.get_bounding_sphere_radius();
     double longest_edge = r/mesh_resolution;
     
     Mesher mesher(cdt);
     if ( !seeds.empty())
          mesher.set_seeds(seeds.begin(),seeds.end());
        
     mesher.set_criteria(Criteria(0.125, longest_edge), true );

     std::cout << "Start  meshing" << std::endl;
 
     mesher.refine_mesh(); // error step_by_step_refine_mesh(); 	
     //-------------------------------------------------------------
     // Remove facets outside can lead to errors if not a simple connected closed polyline.
     //-------------------------------------------------------------
     std::cout << "Done  meshing" << std::endl;

     for(CDT::Face_iterator fit = cdt.faces_begin(); fit != cdt.faces_end(); ++fit)
     {
         if (!fit->is_in_domain())
         {
            cdt.delete_face(fit);
         }
     } 
 


}


void Slice::save(std::string outpath)
{

     if ( cdt.number_of_faces()==0 ) 
     {
        std::cout <<"The resulting mesh has no facet, and will not be saved"<< std::endl;
        return ;     
     }

     std::string extension = outpath.substr(outpath.find_last_of(".")+1);
     std::ofstream out(outpath);
     if ( extension=="off")
     {
       std::ofstream out(outpath);
       CGAL::export_triangulation_2_to_off(out,cdt);
     } 
     else if ( extension=="stl")
     {
        write_STL(outpath);
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
