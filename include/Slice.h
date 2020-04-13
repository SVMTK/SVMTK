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

template <class CDT>
void output_slice_to_medit_(std::ostream& os,
                const CDT& cdt)
{
  // Based on CGAL output_to_medit, but writes all facets to file.
  // i.e. boundary factes and internal facets 
 //TODO ; improve
  typedef typename CDT::Vertex_handle Vertex_handle;
  typedef typename CDT::Face_circulator Face_circulator;
  typedef typename CDT::Face_iterator Face_iterator;
  typedef typename CDT::Vertex_iterator Vertex_iterator;
  typedef typename CDT::Tds Tds;
  typedef typename CDT::Segment Segment;
  typedef typename CDT::Face_circulator Face_circulator;
  typedef typename CDT::Edge_circulator Edge_circulator;
  Tds tds = cdt.tds();

  os << std::setprecision(17);
  os << "MeshVersionFormatted 1\n"
     << "Dimension 2\n";
  //-------------------------------------------------------
  // Vertices
  //-------------------------------------------------------

  os << "Vertices\n" << cdt.number_of_vertices() << '\n';
  boost::unordered_map<Vertex_handle, int> V;
  boost::unordered_map<Vertex_handle, int> E;
  int inum = 1;
  for( Vertex_iterator vit = cdt.vertices_begin();
       vit != cdt.vertices_end();
       ++vit)
  {
    V[vit] = inum++;
    vit->info() = 0; 
    os << *vit <<" "<< 0 <<'\n';
  }
  
  int id=1;

  for(typename CDT::Constraint_iterator cit = cdt.constraints_begin(); cit != cdt.constraints_end(); ++cit)
  {

      for(Vertex_handle vh : cdt.vertices_in_constraint(*cit)) 
      {
        vh->info()=id;
      }
      ++id;
  }
  
  std::set<std::pair<Vertex_handle,Vertex_handle>> set_edges;
  for(Face_iterator ib = cdt.finite_faces_begin(); ib != cdt.finite_faces_end(); ++ib) 
  {
     
     for( int i =0 ; i<3 ; ++i)
     {
         Vertex_handle vh1 = ib->vertex(cdt.ccw(i));
         Vertex_handle vh2 = ib->vertex(cdt.cw(i));
         if (V[vh1] > V[vh2])
         {
            set_edges.insert(std::pair<Vertex_handle,Vertex_handle>(vh1,vh2));
         } 
         else
         {
            set_edges.insert(std::pair<Vertex_handle,Vertex_handle>(vh2,vh1));
         }
     }

   
  }
  os << "Edges\n" 
  << set_edges.size() << '\n';
  for( auto eit : set_edges) 
  {     

         os << V[eit.first] << " " << V[eit.second]  <<" ";
         if  (eit.first->info()==eit.second->info())
         {          
             os << eit.first->info() << std::endl;
         }
         else
         {
             os << 0 << std::endl;
         }
         
     
  }
  os << "Triangles\n" 
  << cdt.number_of_faces()<< '\n';
  for(Face_iterator ib = cdt.finite_faces_begin(); ib != cdt.finite_faces_end(); ++ib) 
  {
     
        os << V[ib->vertex(0)]<< " " <<V[ib->vertex(1)]<<" " <<V[ib->vertex(2)]<< " " <<ib->info() <<"\n";
   
  }





  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End\n";

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
      
       Slice(){}
       ~Slice(){}
 
       Slice(const Polylines_2 &polylines) ;


       void write_STL(const std::string filename);
       void create_mesh(double mesh_resolution);               
       void simplify(double stop_crit); 
       void set_plane(Plane_3 inplane){ this->plane = inplane;}

       void add_subdomains(Domain& domain);
       void save(std::string outpath);

       void add_constraints(Polylines_2 &polylines, bool hole=false); 
       void add_constraints(Slice &slice, bool hole=false);
       void add_constraint(Polyline_2 &polyline, bool hole=false);

       void set_constraints() ;

       void clear_costraints(){constraints.clear();} 

       Polylines_2& get_constraints() { return constraints;}

       std::shared_ptr<Surface> export_3D(); //Forward called

       int  num_constraints() { return constraints.size();}

       // TODO : make functions operational after usage of  polyline

    private:
       Minimum_sphere_2<Kernel> min_sphere;
       Polylines_2 constraints;
       CDT cdt;
       Plane_3 plane;


};


Slice::Slice(const Polylines_2 &polylines) 
{
    min_sphere.add_polylines(polylines);

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
     if (hole)
     {
  
         for ( auto pol =  polylines.begin(); pol != polylines.end(); ++pol)
         {
             Polyline_2 temp = *pol;
             if  ( static_cast<int>(CGAL::orientation_2(temp.begin(),temp.end()))== -1) 
             {
               std::reverse(temp.begin(),temp.end());
             }
             Point_2 c2 = CGAL::centroid(temp.begin(), temp.end(), CGAL::Dimension_tag<0>());
 
             constraints.push_back(temp);

         }
     }
     else 
     {
         constraints.insert( constraints.end(),polylines.begin(),polylines.end());
     }

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
     
     for ( auto  pol : constraints ) 
     {         
               cdt.insert_constraint(pol.begin(), pol.end(),false); 
             
     }
}


void Slice::create_mesh(double mesh_resolution) 
{
     std::cout << "Setting constraints" << std::endl;
     set_constraints(); 
     double r = min_sphere.get_bounding_sphere_radius();
     double longest_edge = r/mesh_resolution;
     
     Mesher mesher(cdt);    
     
     mesher.set_criteria(Criteria(0.125, longest_edge), true );

     std::cout << "Start  meshing "  << std::endl;
 
     mesher.refine_mesh(); 	

     std::cout << "Done  meshing" << std::endl;


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
     else if ( extension=="vtu")
     {
        std::ofstream out(outpath);
        CGAL::write_vtu(out,cdt);
     }
     else if ( extension=="mesh")
     {
         std::ofstream out(outpath);
         output_slice_to_medit_(out,cdt) ;
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
