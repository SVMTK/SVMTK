// Copyright (C) 2018-2020 Lars Magnus Valnes 
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
#ifndef __DOMAIN_H


#define __DOMAIN_H



// LOCAL
#include "Surface.h" 
#include "Slice.h" 
#include "SubdomainMap.h" 
#include "Polyhedral_vector_to_labeled_function_wrapper.h"
//#include "Labeled_mesh_domain_with_features.h"
// STD 
#include <list>
#include <fstream>
#include <memory>
#include <set> 
#include <map>
#include <fstream>
//CGAL MESH_3
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // CGALSURFACE  CGALSURFACE ? 
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_3/Detect_polylines_in_polyhedra.h>
#include <CGAL/Mesh_3/polylines_to_protect.h>
#include <CGAL/Mesh_3/Mesher_3.h>

#include <CGAL/Mesh_3/C3T3_helpers.h>

#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// CGAL BOUNDING SPHERE
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/internal/Mesh_3/Boundary_of_subdomain_of_complex_3_in_triangulation_3_to_off.h>
#include <CGAL/Triangulation_3.h>
// CURRENT PLAN :
// TODO: Consider void remove_isolated_vertices() as stand alone function or class function
// TODO: Clean up header files

template<class C3T3, class PointContainer, class FaceContainer>
void facets_in_complex_3_to_triangle_soup_(const C3T3& c3t3,
                                          const typename C3T3::Surface_patch_index sf_index,
                                          PointContainer& points,
                                          FaceContainer& faces,
                                          const bool normals_point_outside_of_the_subdomain = false,
                                          const bool export_all_facets = false)
                                       
{
  typedef typename PointContainer::value_type                            Point_3;
  typedef typename FaceContainer::value_type                             Face;

  typedef typename C3T3::Triangulation                                   Tr;

  typedef typename Tr::Vertex_handle                                     Vertex_handle;
  typedef typename Tr::Cell_handle                                       Cell_handle;
  typedef typename Tr::Weighted_point                                    Weighted_point;

  typedef typename C3T3::Facets_in_complex_iterator                      Ficit;

  typedef CGAL::Hash_handles_with_or_without_timestamps                  Hash_fct;
  typedef boost::unordered_map<Vertex_handle, std::size_t, Hash_fct>     VHmap;

  typedef typename C3T3::size_type                                       size_type;

  size_type nf = c3t3.number_of_facets_in_complex();
  faces.reserve(faces.size() + nf);
  points.reserve(points.size() + nf/2); // approximating Euler

  VHmap vh_to_ids;
  std::size_t inum = 0;
  int df =0;
  for(Ficit fit = c3t3.facets_in_complex_begin(sf_index),
            end = c3t3.facets_in_complex_end(sf_index); fit != end; ++fit)
  {
    Cell_handle c = fit->first;
    int s = fit->second;
    Face f;
    CGAL::Mesh_3::internal::resize(f, 3);

    typename C3T3::Subdomain_index cell_sdi = c3t3.subdomain_index(c);
    typename C3T3::Subdomain_index opp_sdi = c3t3.subdomain_index(c->neighbor(s));

    //if(cell_sdi != sf_index.second or opp_sdi != sf_index.first)
    //   continue;
    
    for(std::size_t i=1; i<4; ++i)
    {
      typename VHmap::iterator map_entry;
      bool is_new;
      Vertex_handle v = c->vertex((s+i)&3);
      CGAL_assertion(v != Vertex_handle() && !c3t3.triangulation().is_infinite(v));

      boost::tie(map_entry, is_new) = vh_to_ids.insert(std::make_pair(v, inum));
      if(is_new)
      {
        const Weighted_point& p = c3t3.triangulation().point(c, (s+i)&3);
        const Point_3 bp = Point_3(CGAL::to_double(p.x()),
                                   CGAL::to_double(p.y()),
                                   CGAL::to_double(p.z()));
        points.push_back(bp);
        ++inum;
      }

      f[i-1] = map_entry->second;
    }
    if( sf_index.first < sf_index.second)
        std::swap(f[0], f[1]);
    faces.push_back(f);
  }
  
}


template <class C3T3,
          class Vertex_index_property_map,
          class Facet_index_property_map,
          class Facet_index_property_map_twice,
          class Cell_index_property_map>
void output_to_medit_(std::ostream& os,
                const C3T3& c3t3,
                const Vertex_index_property_map& vertex_pmap,
                const Facet_index_property_map& facet_pmap,
                const Cell_index_property_map& cell_pmap,
                const Facet_index_property_map_twice& facet_twice_pmap = Facet_index_property_map_twice(),
                const bool print_each_facet_twice = false,
                const bool save_edges = false)
{
  // Based on CGAL output_to_medit, but writes all facets to file.
  // i.e. boundary factes and internal facets 
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Tr::Cell_circulator Cell_circulator;

  const Tr& tr = c3t3.triangulation();

  os << std::setprecision(17);

  os << "MeshVersionFormatted 1\n"
     << "Dimension 3\n";

  os << "Vertices\n" << tr.number_of_vertices() << '\n';

  boost::unordered_map<Vertex_handle, int> V;
  int inum = 1;
  for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
       vit != tr.finite_vertices_end();
       ++vit)
  {
    V[vit] = inum++;
    Weighted_point p = tr.point(vit);
    os << CGAL::to_double(p.x()) << ' '
       << CGAL::to_double(p.y()) << ' '
       << CGAL::to_double(p.z()) << ' '
       << get(vertex_pmap, vit)
       << '\n';
  }

  if ( save_edges) 
  {
     bool flag=false;
     int number_of_edges=0;
     // TODO: improve
     for ( auto eit = tr.edges_begin() ; eit != tr.edges_end(); ++eit) 
     {       
         flag=false; //TODO: use break or continue??


         Cell_circulator ccir = tr.incident_cells(*eit);
         Cell_circulator cdone = ccir;
         do 
         {
             if ( c3t3.is_in_complex(ccir) )
             { 
                flag=true; 

             }
         *ccir++;
         }while(ccir!=cdone and flag==false);          

         if( flag ) 
         {
            number_of_edges++;
         }
     }
     
     os << "Edges\n" 
     << number_of_edges << '\n';
     for( auto eit = tr.finite_edges_begin(); eit != tr.finite_edges_end(); ++eit) 
     {
         flag=false;
         Cell_circulator ccir = tr.incident_cells(*eit);
         Cell_circulator cdone = ccir;
         do 
         {
             if ( c3t3.is_in_complex(ccir) )
             { 
                flag=true; 

             }
         *ccir++;
         }while(ccir!=cdone and flag==false);          

         if( flag ) 
         {         
            Vertex_handle vh1 = eit->first->vertex(eit->second);
            Vertex_handle vh2 = eit->first->vertex(eit->third);
            os << V[vh1] << " " << V[vh2]  <<" " << c3t3.curve_index(*eit) << std::endl;
         }
     }

  }

  //-------------------------------------------------------
  // Facets
  //-------------------------------------------------------
  int number_of_triangles=0;
  for( auto fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)  
  {
      if (c3t3.is_in_complex(fit->first) or  c3t3.is_in_complex(fit->first->neighbor(fit->second) ) ) 
      number_of_triangles++;
  }

  os << "Triangles\n" 
     << number_of_triangles << '\n';

  for( auto fit = tr.finite_facets_begin(); fit != tr.finite_facets_end(); ++fit)
  {
    typename C3T3::Facet f = (*fit);
    
      // Apply priority among subdomains, to get consistent facet orientation per subdomain-pair interface.

      // NOTE: We mirror a facet when needed to make it consistent with No_patch_facet_pmap_first/second.
    if (f.first->subdomain_index() > f.first->neighbor(f.second)->subdomain_index())
    {
        f = tr.mirror_facet(f);
    }
    
    // Get facet vertices in CCW order.
    Vertex_handle vh1 = f.first->vertex((f.second + 1) % 4);
    Vertex_handle vh2 = f.first->vertex((f.second + 2) % 4);
    Vertex_handle vh3 = f.first->vertex((f.second + 3) % 4);
    
    // Facet orientation also depends on parity.
    if (f.second % 2 != 0)
      std::swap(vh2, vh3);
    
    if (c3t3.is_in_complex(fit->first) or  c3t3.is_in_complex(fit->first->neighbor(fit->second) ) )  {

    os << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3] << ' '; 
    if (get(facet_pmap, *fit)<0)
    {
      os << 0 << '\n'; // adjust -1 -> 0
    }
    else 
    {
      os << get(facet_pmap, *fit) << '\n'; // adjust -1 -> 0

    } 

    }


  }

  //-------------------------------------------------------
  // Tetrahedra
  //-------------------------------------------------------
  os << "Tetrahedra\n"
     << c3t3.number_of_cells_in_complex() << '\n';

  for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit )
  {
    for (int i=0; i<4; i++)
      os << V[cit->vertex(i)] << ' ';

    os << get(cell_pmap, cit) << '\n';
  }

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End\n";

}
template<typename Kernel>
struct Minimum_sphere
{
   
           typedef typename CGAL::Min_sphere_of_spheres_d_traits_3<Kernel, typename Kernel::FT> MinSphereTraits;
           typedef typename CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
           typedef typename MinSphereTraits::Sphere Sphere;

           template< typename MeshPolyhedron_3>  
           void add_polyhedron(const MeshPolyhedron_3 &polyhedron)
           {
                for (typename MeshPolyhedron_3::Vertex_const_iterator it=polyhedron.vertices_begin();it != polyhedron.vertices_end(); ++it)
                {
                    S.push_back(Sphere(it->point(), 0.0));
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




class Surface;
class Slice; 
 

class Domain {
    public :

        //------- Kernel ------------
        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

        typedef Kernel::Point_3 Point_3;
      
        typedef CGAL::Mesh_polyhedron_3<Kernel>::type Polyhedron; 
        
        typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel, Polyhedron> Polyhedral_mesh_domain_3; 
        typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, Kernel  > Function_wrapper; 

        typedef Function_wrapper::Function_vector Function_vector; 
        typedef CGAL::Labeled_mesh_domain_3<Kernel> Labeled_Mesh_Domain;
        typedef CGAL::Mesh_domain_with_polyline_features_3<Labeled_Mesh_Domain> Mesh_domain; 


        // CGAL::Parallel_if_available_tag
        #ifdef CGAL_CONCURRENT_MESH_3
               typedef CGAL::Parallel_tag Concurrency_tag;
        #else
               typedef CGAL::Sequential_tag Concurrency_tag;
        #endif
        typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

        typedef int Curve_index; // FIXME  ?
        typedef int Corner_index; //FIXME ?

        // experimental adjustsments to class ; Corner_index and Curve_index
        typedef CGAL::Mesh_complex_3_in_triangulation_3< Tr, Corner_index, Curve_index> C3t3;


        typedef C3t3::Subdomain_index Subdomain_index;
        typedef C3t3::Cell_handle Cell_handle;
        typedef C3t3::Vertex_handle Vertex_handle;
        typedef C3t3::Surface_patch_index Surface_patch_index;
        typedef C3t3::Cells_in_complex_iterator Cell_iterator;
        typedef C3t3::Facets_in_complex_iterator Facet_iterator;
        
        typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
        typedef CGAL::Triple<Cell_handle, int, int> Edge; 


        // ------ Auxiliary, also  experimental----------------
        typedef Tr::Finite_vertices_iterator Finite_vertices_iterator;
        typedef Tr::Locate_type Locate_type;
        typedef Tr::Weighted_point Weighted_point;
        typedef Tr::Facet Facet;
        typedef std::vector<Point_3>  Polyline_3;
        typedef std::vector<Polyline_3> Polylines;
        typedef std::map<std::string, double> Parameters;     
        typedef std::vector<std::size_t>  Face;


     

        // --------Constructors---------------------
        Domain(Surface& surface);
        Domain(std::vector<Surface> surfaces);
        Domain(std::vector<Surface> surfaces, AbstractMap& map);
        // -----------------------------------------
        ~Domain() {}

        
        // --------CGAL mesh criteria---------
        void set_parameters(std::map<std::string, double> &new_parameters); // py::dict? 
        void set_parameter(std::string key, double value);

        // ------- Mesh Creation and Refinement
        void create_mesh(const double mesh_resolution );
        void create_mesh();
        void refine_mesh(const double mesh_resolution );
        void refine_mesh();
        // --------CGAL mesh critera settings
        void default_parameters() {
            parameters["mesh_resolution"]=64.0;
            parameters["facet_angle"]    = 25.0;
            parameters["facet_size"]     = 0.1;
            parameters["facet_distance"] =  0.1;
            parameters["cell_radius_edge_ratio"] = 3.0;
            parameters["cell_size"] = 0.1;
            parameters["edge_size"] = 0.1;
        }
        // ----  IO Funtions
        void save(std::string OutPath, bool save_1Dfeatures); 



        void add_surface_points(Surface &surface);
        void add_sharp_border_edges(Polyhedron& polyhedron, double threshold);
        void add_sharp_border_edges(Surface& surface, double threshold=60); 
 

        Polylines& get_features() {return features; }  
        Polylines& get_borders() {return borders; }
        void add_feature(Polyline_3 polyline) { features.push_back(polyline);}  

       
        void clear_borders(){borders.clear();} // rename

        void set_borders(){domain_ptr.get()->add_features(get_borders().begin(), get_borders().end());} // -> TODO: integrate or REOMVE 
        void add_borders(Polyline_3 polyline) { borders.push_back(polyline);} 
      

   
        std::set<std::pair<int,int>> number_of_patches();
        void remove_subdomain(std::vector<int> tags);
        void remove_subdomain(int tag);              // tags 

        void set_features();
        void set_features(Polylines& polylines){ set_features( polylines.begin() , polylines.end() );} 

        template< typename InputIterator> // Wrapping 
        void set_features(InputIterator begin, InputIterator end){domain_ptr->add_features(begin, end);} 

        std::shared_ptr<Slice> mesh_slice(double x1,double x2, double x3 ,double x4);
        std::shared_ptr<Surface> get_boundary(int tag); // # subdomainindex


        // -----Mesh optimization ->Direct wrapping of CGAL funcions
        void lloyd(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true);// default values
        void odt(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true);  // default values
        void exude( double time_limit = 0, double sliver_bound = 0 ){ CGAL::exude_mesh_3(c3t3, sliver_bound=sliver_bound, time_limit=time_limit);}  // default values
        void perturb( double time_limit=0, double sliver_bound=0){CGAL::perturb_mesh_3 ( c3t3, *domain_ptr.get(), time_limit=time_limit, sliver_bound=sliver_bound) ;} // default values



        void rebind_facets();
        int number_of_cells(){return c3t3.number_of_cells();}
        std::set<int>  number_of_subdomains();
        C3t3& get_mesh(){return c3t3;}
 

    private :
        std::unique_ptr<Mesh_domain> domain_ptr;
        Minimum_sphere<Kernel> min_sphere; 
        Parameters parameters; 
        C3t3 c3t3;
        Polylines borders; 
        Polylines features;

};



template<typename C3T3> // TODO: mesh function 
int remove_isolated_vertices(C3T3& c3t3, bool remove_domain=false)
{ 
  //FIXME : CLEANUP
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;

  std::map<Vertex_handle, bool> vertex_map;

  for( Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin();vit != c3t3.triangulation().finite_vertices_end();++vit)
  { 
        vertex_map[vit] = false ;  
  }

  for(Cell_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
    for (std::size_t i = 0; i < 4; ++i) // 
    {
        vertex_map[cit->vertex(i)] = true;
    }
  }

  int before = c3t3.triangulation().number_of_vertices() ;

  for (typename std::map<Vertex_handle, bool>::const_iterator it = vertex_map.begin();it != vertex_map.end(); ++it) // check post or pre increment
  {
    if (!it->second) 
    {
       c3t3.triangulation().remove(it->first);
    }

  }

  int after = c3t3.triangulation().number_of_vertices() ; 

  std::cout<<"Number of isolated vertices removed: "<< before - after << std::endl;
  if ((before - after) > 10 and !remove_domain)   
  {
       // FIXME : Better 
       std::cout<<"The number of removed vertices are substantial."<< std::endl;  
       std::cout<<"This can cause the mesh to be horrid."<< std::endl;  
       std::cout<<"Try isotropic remeshing or increase the mesh resolution."<< std::endl;  
  }

  return (before - after);


}


void Slice::add_subdomains(Domain& domain)
{
   typedef  Domain::C3t3  C3t3;
   typedef  Domain::Point_3 Point_3; 
   typedef  Domain::Tr Tr; 
   typedef  Domain::Cell_handle Cell_handle;
   typedef  Domain::Weighted_point Weighted_point;
   typedef  Domain::Surface_patch_index Surface_patch_index;
   typedef  Domain::Tr::Locate_type  Locate_type;
   typedef  Domain::Facet Facet;
   typedef  std::pair<int,int>  Pid; // patch id 
   typedef  std::map<Pid,int> Pid_map;
  

   int j, fi,fn,ci,cn;
   Cell_handle cin,cit; 
   Pid_map pid_map_;
   Pid spp;
   int index_counter=1;
   Domain::C3t3 c3t3 = domain.get_mesh(); 
     //auto patches = domain.number_of_patches() ;

   Tr tr= c3t3.triangulation();

   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
      Point_2 c2 =  CGAL::centroid(fit->vertex(0)->point() ,fit->vertex(1)->point(),fit->vertex(2)->point());  
      Point_3 c3 = this->plane.to_3d(c2);
      Weighted_point wp = Weighted_point(c3,1.0) ;
      cit = tr.inexact_locate(wp);
      ci = static_cast<int>(c3t3.subdomain_index(cit));
     
      fit->info() = ci;
   }
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {
          fi = fit->info(); 
          for( int i =0 ;i<3;++i)
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
                  std::cout <<  spp.first << "\t" << spp.second << index_counter<< std::endl;
                }
                this->edges[ei] = pid_map_[spp];    
             }
             else
             {
                this->edges[ei] =  0;
         
             }
           
         }   
        
    }

   j=0;
   Locate_type loc;
   
   int li,lj;
   for(Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
   {


      for( int i =0 ;i<3;++i)
      {
         Vertex_handle vh1 = fit->vertex(cdt.ccw(i));
         Vertex_handle vh2 = fit->vertex(cdt.cw(i));
         
         Point_3 pn3 =  this->plane.to_3d(vh1->point());
         Weighted_point wpn = Weighted_point(pn3,1.0) ;
         Cell_handle cin = tr.locate(wpn,loc,li,lj);
         Domain::Surface_patch_index spi = c3t3.surface_patch_index(cin,li);
         //std::cout<< spi.first <<"\t" << spi.second << std::endl;
         //
         //
           
       }
    }

/*
             std::cout<< j<< std::endl;
             if (cn!=ci) 
             { 
                if (ci>cn)
                    spp={ci , cn};
                else 
                    spp={cn , ci};

                std::pair<Pid_map::iterator, bool> is_insert_successful = pid_map_.insert(std::make_pair(spp,index_counter));
                if(is_insert_successful.second)
                    ++index_counter;

                this->edges[eit] = pid_map_[spp] ;
             }
             else  
             {
                this->edges[eit] = 0;              
             }
 
          

      }
 
   }*/ 


} 
// --------------------------------
//  Constructor for Surfaces  
//  TODO: smart pointer error : 
//  python3: ../nptl/pthread_mutex_lock.c:352: __pthread_mutex_lock_full: Assertion `INTERNAL_SYSCALL_ERRNO (e, __err) != ESRCH || !robust' failed
// --------------------------------
Domain::Domain( std::vector<Surface> surfaces )
{
    Function_vector v;;
    for(std::vector<Surface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       Polyhedron *polyhedron = new Polyhedron();
       sit->get_polyhedron(*polyhedron);
       min_sphere.add_polyhedron(*polyhedron);

       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(*polyhedron);
       v.push_back(polyhedral_domain);
    }
    Function_wrapper wrapper(v);
    domain_ptr = std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));
    default_parameters();
}

Domain::Domain( std::vector<Surface> surfaces , AbstractMap& map )
{
    Function_vector v;
    
    for(std::vector<Surface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       Polyhedron *polyhedron = new Polyhedron();
       sit->get_polyhedron(*polyhedron);
       min_sphere.add_polyhedron(*polyhedron);
       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(*polyhedron);
       v.push_back(polyhedral_domain);
      
    }
    Function_wrapper wrapper(v,map);

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));
    default_parameters();
}

Domain::Domain(Surface &surface) 
{
    Function_vector v;
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);

    min_sphere.add_polyhedron(polyhedron);

    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);

    v.push_back(polyhedral_domain);
    Function_wrapper wrapper(v);

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox())); 
    default_parameters();

}

// --------------------------------
//  
// --------------------------------

std::shared_ptr<Surface> Domain::get_boundary(int tag=0)
{
   std::vector<Face> faces;
   Polyline_3 points;
    
   Subdomain_index useless = Subdomain_index(tag);
   facets_in_complex_3_to_triangle_soup(c3t3, useless, points, faces, true, false);  

   std::shared_ptr<Surface> surf(new  Surface(points, faces)) ;
   return surf;
}


// ---------- CGAL mesh criteria --------------------
void Domain::set_parameters(std::map<std::string, double> &new_parameters)
{
    for (std::map<std::string, double>::iterator pit= new_parameters.begin(); pit!=new_parameters.end(); ++pit )
    {
        parameters[pit->first] = static_cast<double>(pit->second);
    }
}
void Domain::set_parameter(std::string key , double value )
{
   parameters[key] = value;
}

void Domain::create_mesh()
{
    std::cout << "begin_meshing" << std::endl;
 
    Mesh_criteria criteria(CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );

    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);

    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();


}

void Domain::create_mesh(const double mesh_resolution )
{

    double r = min_sphere.get_bounding_sphere_radius(); 
    const double cell_size = r/mesh_resolution;


    std::cout << "Cell size: " << cell_size << std::endl;
    Mesh_criteria criteria(CGAL::parameters::edge_size = cell_size,
                                       CGAL::parameters::facet_angle = 30.0,
                                       CGAL::parameters::facet_size = cell_size,
                                       CGAL::parameters::facet_distance = cell_size/10.0, 
                                       CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                       CGAL::parameters::cell_size = cell_size);

    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);


    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();
  
}

void Domain::save(std::string OutPath,bool save_1Dfeatures)
{
    
    std::ofstream  medit_file(OutPath);


    typedef CGAL::Mesh_3::Medit_pmap_generator<C3t3,false,false> Generator;
    typedef typename Generator::Cell_pmap Cell_pmap;
    typedef typename Generator::Facet_pmap Facet_pmap;
    typedef typename Generator::Facet_pmap_twice Facet_pmap_twice;
    typedef typename Generator::Vertex_pmap Vertex_pmap;

    Cell_pmap cell_pmap(c3t3);
    Facet_pmap facet_pmap(c3t3, cell_pmap);
    Facet_pmap_twice  facet_twice_pmap(c3t3,cell_pmap);
    Vertex_pmap vertex_pmap(c3t3, cell_pmap,facet_pmap);

    output_to_medit_(medit_file,c3t3, vertex_pmap, facet_pmap, cell_pmap, facet_twice_pmap , false, save_1Dfeatures) ;
    medit_file.close();
}


void Domain::refine_mesh()
{
    Mesh_criteria criteria(CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );

   refine_mesh_3(c3t3, *domain_ptr.get(), criteria,CGAL::parameters::no_reset_c3t3());

    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();

}
void Domain::refine_mesh(const double mesh_resolution )
{
    double r = min_sphere.get_bounding_sphere_radius(); 
    const double cell_size = r/mesh_resolution;

    Mesh_criteria criteria(CGAL::parameters::edge_size = cell_size,
                                       CGAL::parameters::facet_angle = 30.0,
                                       CGAL::parameters::facet_size = cell_size,
                                       CGAL::parameters::facet_distance = cell_size/10.0, 
                                       CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                       CGAL::parameters::cell_size = cell_size);

   refine_mesh_3(c3t3, *domain_ptr.get(), criteria,CGAL::parameters::no_reset_c3t3());

   if ( remove_isolated_vertices(c3t3) > 0)
       c3t3.rescan_after_load_of_triangulation();

}
void Domain::add_surface_points(Surface &surface) 
{ 

   typedef Tr::Bare_point Bare_point;
   typedef std::vector<std::pair<Bare_point, int> > Initial_points_vector;
   typedef Surface::Index SM_Index;
   std::vector<Point_3> points;  
   std::vector<SM_Index> indicies; 
   Initial_points_vector points_vector;
   points = surface.get_points() ; 
   indicies = surface.get_vertices(); 
   for(int i = 0; i<points.size(); i++)
   {
      points_vector.push_back(std::make_pair<Bare_point,int>(Bare_point(points[i]), indicies[i] )    ) ;
   }
   c3t3.insert_surface_points(points_vector.begin(), points_vector.end());
}


void Domain::remove_subdomain(int tag) //rename
{ 
   std::vector<int> temp;
   temp.push_back(tag);
   remove_subdomain(temp);
}




void Domain::remove_subdomain(std::vector<int> tags) //rename
{

   // TODO : CLEAN UP 
  int before = c3t3.number_of_cells();

  std::map<Cell_handle,int> cell_map;
  std::vector<std::tuple<Cell_handle,int,int>> rebind;

  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  { 
     cell_map[cit]=static_cast<int>( c3t3.subdomain_index(cit) );
  }

  for( auto j = tags.begin(); j!=tags.end(); ++j)
  {
    Subdomain_index temp(*j) ; 
    for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
    {
      for (std::size_t i = 0; i < 4; i++)
      { 
         if(std::find(tags.begin(), tags.end(), static_cast<int>(c3t3.subdomain_index(cit))  ) == tags.end() ) // not a cell to be removed               
         {
           if(c3t3.subdomain_index(cit->neighbor(i))==temp)
           {
               rebind.push_back(std::make_tuple(cit, i,*j)); 
           }
         }
      }
    }
  }
     
  for(std::vector<int>::iterator j = tags.begin(); j!=tags.end(); ++j)
  {
      for(auto cit =  c3t3.cells_in_complex_begin(Subdomain_index(*j));cit !=  c3t3.cells_in_complex_end(); ++cit)
      {  
          c3t3.remove_from_complex(cit); 
      }     
  }


  c3t3.rescan_after_load_of_triangulation(); 
  remove_isolated_vertices(c3t3,true);

  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
      for(int i=0; i<4; ++i)
      {
        Cell_handle cn = cit;
        Cell_handle cm = cit->neighbor(i);
        Subdomain_index ci = c3t3.subdomain_index(cn);
        Subdomain_index cj = c3t3.subdomain_index(cm);
       
        if( ci!=cj )
        {   
              if ( ci > cj ) 
              {
              c3t3.remove_from_complex(cit,i);
              c3t3.add_to_complex(cit, i, Surface_patch_index(ci,cj) ); 
              }   
              else 
              {
              c3t3.remove_from_complex(cit,i);
              c3t3.add_to_complex(cit, i, Surface_patch_index(cj,ci) ); 
              }
        }
      }
  }

  for ( auto cit = rebind.begin(); cit!=rebind.end(); ++cit)
  {
      Cell_handle cn = std::get<0>(*cit);
      int s =std::get<1>(*cit);
      Cell_handle cm=  cn->neighbor(s);
      Subdomain_index ck = Subdomain_index(std::get<2>(*cit));            
      c3t3.remove_from_complex(cn,s);
      c3t3.add_to_complex(cn,s,Surface_patch_index(ck,0) ); //change to 
  }

  int after = c3t3.number_of_cells();

  std::cout << "Number of removed subdomain cells : " << (before -after) << std::endl;
  c3t3.rescan_after_load_of_triangulation(); //why??
}



void Domain::add_sharp_border_edges(Polyhedron& polyhedron, double threshold)
{ 
 Polylines polylinput; 
 typedef boost::property_map<Polyhedron, CGAL::edge_is_feature_t>::type EIF_map;
 EIF_map eif = get(CGAL::edge_is_feature, polyhedron);

 CGAL::Polygon_mesh_processing::detect_sharp_edges(polyhedron,threshold, eif); 
 for( Polyhedron::Edge_iterator he = polyhedron.edges_begin(); he != polyhedron.edges_end() ; ++he)
 {
   if(he->is_feature_edge() ) 
   {
       Polyline_3 polyline;
       polyline.push_back(he->vertex()->point());
       polyline.push_back(he->opposite()->vertex()->point());     
       borders.push_back(polyline);
    }    
 }   
}

void Domain::add_sharp_border_edges(Surface& surface, double threshold) 
{ 
  Polyhedron polyhedron;
  surface.get_polyhedron(polyhedron);
  add_sharp_border_edges(polyhedron, threshold);
}

std::set<int>  Domain::number_of_subdomains()
{
   std::set<int> sd_indices;
   for(Cell_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
   {
        sd_indices.insert(static_cast<int>(c3t3.subdomain_index(cit)) );
   }
   return sd_indices;
}
std::set<std::pair<int,int>>  Domain::number_of_patches()
{
   std::set<std::pair<int,int>> sf_indices;
   for(Cell_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
   {
        for( int i =0 ; i<4 ; ++i)
        {
        Surface_patch_index spi = c3t3.surface_patch_index(cit,i) ;
        sf_indices.insert(std::pair<int,int>(static_cast<int>(spi.first) , static_cast<int>(spi.second) ));
        }
   }
   return sf_indices;
}

//() ? ?
std::shared_ptr<Slice> Domain::mesh_slice(double x1,double x2, double x3 ,double x4) 
{ 
   Kernel::Plane_3 plane = Kernel::Plane_3(x1, x2, x3, x4);

   std::shared_ptr<Slice> slice(new  Slice(plane)) ;
   //slice->set_plane(plane);

   std::vector<Point_3> points;
   std::vector<Face> faces;

   for ( auto i : number_of_patches() ) 
   {
      points.clear(); 
      faces.clear();
      if (i.first != i.second)
      { 
      facets_in_complex_3_to_triangle_soup_(c3t3,Surface_patch_index(i.first,i.second),points,faces);
      std::cout << "out \t"<< faces.size() << "\t" << i.first << "\t" << i.second << std::endl;
   
      std::shared_ptr<Surface> surf(new Surface(points,faces)); 

      std::shared_ptr<Slice> temp = surf->mesh_slice(x1,x2,x3,x4);     // nerror here   

      slice->add_constraints(*temp.get());  // seg fault ??
    
      surf->save("surf"+std::to_string(i.first)+std::to_string(i.second)+".off");
      }
   }

   return slice;
}
/*
void Domain::load_medit(std::string filename )
{

   // TODO: without tokenizer ?? 
   typedef boost::tokenizer<boost::char_separator<char> > tokenizer; 

   const boost::char_separator<char> sep(" ");
   Tr& tr = c3t3.triangulation();
   //assert(tr.num_vertices()==0);

   std::ifstream in(filename,std::ios_base::in);
   std::map<Vertex_handle,int> V; 
   std::map<Cell_handle, std::vector<int> > F; 


   std::string line;
   int linenum = 0
   int num_vertices = 0
   double x1,x2,x3,x4; 
   int v1,v2,v3,v4;

   int point=0;
   while (std::getline(in, line))
   {
       linenum++; 
       if (linenum==1)
          assert(line.compare("MeshVersionFormatted 1") );

       if (linenum==2) 
          assert(line.compare("Dimension 2") );

       if(linenum==3) 
           assert( line.comapre("Vertices") );

       if(linenum==4) 
          num_vertices = std::stoi(line);


       if(linenum==5) 
       {
          line >> x1 >> x2 >> x3 >> x4 ; 
          Weighted_point p(x1,x2,x3)        
          V[tr.insert(p1)]=point++; 
       }
      

       if ()
       {set_vertices();
  Cell_handle ch = tr.finite_cells_begin();
    c3t3.add_to_complex(ch,subdomain_index);
        
       }


       if(linenum==(5+num_vertices) 
         assert(line.compare("Edges") ) ;  



       std::istringstream iss(line);




}*/


void Domain::set_features(){set_features(features.begin() , features.end() );}


void Domain::lloyd(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze )
{CGAL::lloyd_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 

void Domain::odt(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze) 
{CGAL::odt_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 


#endif
