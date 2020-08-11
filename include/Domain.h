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
#include "Polyhedral_vector_to_labeled_function_wrapper.h"
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
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Mesh_3/C3T3_helpers.h>

#include <CGAL/refine_mesh_3.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/Triangulation_3.h>


template<class C3T3, class PointContainer, class FaceContainer>
void facets_in_complex_3_to_triangle_soup_(const C3T3& c3t3,
                                          const typename C3T3::Surface_patch_index sf_index,
                                          PointContainer& points,
                                          FaceContainer& faces,
                                          const bool normals_point_outside_of_the_subdomain = false,
                                          const bool export_all_facets = false)
                                       
{
  typedef typename PointContainer::value_type         Point_3;
  typedef typename FaceContainer::value_type          Face;
  typedef typename C3T3::Triangulation                Tr;
  typedef typename C3T3::Facets_in_complex_iterator   Ficit;
  typedef typename C3T3::size_type                    size_type;

  typedef typename Tr::Vertex_handle                  Vertex_handle;
  typedef typename Tr::Cell_handle                    Cell_handle;
  typedef typename Tr::Weighted_point                 Weighted_point;

  typedef CGAL::Hash_handles_with_or_without_timestamps                  Hash_fct;
  typedef boost::unordered_map<Vertex_handle, std::size_t, Hash_fct>     VHmap;

                                      

  size_type nf = c3t3.number_of_facets_in_complex();
  faces.reserve(faces.size() + nf);
  points.reserve(points.size() + nf/2); 
  VHmap vh_to_ids;
  std::size_t inum = 0;
  for(Ficit fit = c3t3.facets_in_complex_begin(sf_index),
            end = c3t3.facets_in_complex_end(sf_index); fit != end; ++fit)
  {
    Cell_handle c = fit->first;
    int s = fit->second;
    Face f;
    CGAL::Mesh_3::internal::resize(f, 3);


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

/** Based on CGAL output_to_medit, but writes all facets to file.
 *  i.e. boundary factes and internal facets 
*/
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

  typedef typename C3T3::Triangulation Tr;
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
     for ( auto eit = tr.edges_begin() ; eit != tr.edges_end(); ++eit) 
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
    if (f.first->subdomain_index() > f.first->neighbor(f.second)->subdomain_index())
    {
        f = tr.mirror_facet(f);
    }

    Vertex_handle vh1 = f.first->vertex((f.second + 1) % 4);
    Vertex_handle vh2 = f.first->vertex((f.second + 2) % 4);
    Vertex_handle vh3 = f.first->vertex((f.second + 3) % 4);
    
    if (f.second % 2 != 0)
      std::swap(vh2, vh3);
    
    if (c3t3.is_in_complex(fit->first) or  c3t3.is_in_complex(fit->first->neighbor(fit->second) ) )  {

    os << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3] << ' '; 
    if (get(facet_pmap, *fit)<0)
    {
      os << 0 << '\n'; 
    }
    else 
    {
      os << get(facet_pmap, *fit) << '\n'; 

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


template<typename C3T3> 
int remove_isolated_vertices(C3T3& c3t3, bool remove_domain=false)
{ 
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
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

  for (typename std::map<Vertex_handle, bool>::const_iterator it = vertex_map.begin();it != vertex_map.end(); ++it)   {
    if (!it->second) 
    {
       c3t3.triangulation().remove(it->first);
    }

  }

  int after = c3t3.triangulation().number_of_vertices() ; 

  std::cout<<"Number of isolated vertices removed: "<< before - after << std::endl;
  if ((before - after) > 10 and !remove_domain)   
  {
       std::cout<<"There were a substantial number of isolated vertices, and the user should inspect the mesh."<< std::endl;  
       std::cout<<"Try: 1.isotropic remeshing, 2.increase the mesh resolution or \n 3. specific mesh parameters to decrease number of isolated vertices."<< std::endl;  

  }

  return (before - after);


}

class Domain {
    public :

        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

        typedef Kernel::Point_3 Point_3;
      
        typedef CGAL::Mesh_polyhedron_3<Kernel>::type Polyhedron; 
        
        typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel, Polyhedron> Polyhedral_mesh_domain_3; 
        typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, Kernel  > Function_wrapper; 

        typedef Function_wrapper::Function_vector Function_vector; 
        typedef CGAL::Labeled_mesh_domain_3<Kernel> Labeled_Mesh_Domain;
        typedef CGAL::Mesh_domain_with_polyline_features_3<Labeled_Mesh_Domain> Mesh_domain; 

        typedef CGAL::Sequential_tag Concurrency_tag;

        typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

        typedef int Curve_index; 
        typedef int Corner_index; 

        typedef CGAL::Mesh_complex_3_in_triangulation_3< Tr, Corner_index, Curve_index> C3t3;


        typedef C3t3::Subdomain_index Subdomain_index;
        typedef C3t3::Cell_handle Cell_handle;
        typedef C3t3::Vertex_handle Vertex_handle;
        typedef C3t3::Surface_patch_index Surface_patch_index;
        typedef C3t3::Cells_in_complex_iterator Cell_iterator;
        typedef C3t3::Facets_in_complex_iterator Facet_iterator;
        
        typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
        typedef CGAL::Triple<Cell_handle, int, int> Edge; 

        typedef Tr::Finite_vertices_iterator Finite_vertices_iterator;
        typedef Tr::Locate_type Locate_type;
        typedef Tr::Weighted_point Weighted_point;
        typedef Tr::Facet Facet;

        typedef std::vector<Point_3>  Polyline_3;
        typedef std::vector<Polyline_3> Polylines;
        typedef std::map<std::string, double> Parameters;     
        typedef std::vector<std::size_t>  Face;

        template<typename Surface>
        Domain(Surface& surface);
        template<typename Surface>
        Domain(std::vector<Surface> surfaces);
        template<typename Surface>
        Domain(std::vector<Surface> surfaces, AbstractMap& map);

        ~Domain() { for( auto vit : this->v){delete vit;}v.clear();}        

        void create_mesh(const double mesh_resolution );
        void create_mesh(double cell_size, double facet_size,double facet_angle,  double facet_distance,double cell_radius_edge_ratio);

        void save(std::string OutPath, bool save_1Dfeatures); 

        double get_bounding_sphere_radius(){ return min_sphere.get_bounding_sphere_radius(); }

        void add_sharp_border_edges(Polyhedron& polyhedron, double threshold);
        template<typename Surface>
        void add_sharp_border_edges(Surface& surface, double threshold=60);
        //template<typename Surface> 
        //void add_surface_points(Surface &surface); // FIXME operator ?? 

        void set_borders();
        void set_features();

        void clear_borders(){this->borders.clear();}
        void clear_features(){this->features.clear();}
        template<typename Surface>
        std::vector<std::shared_ptr<Surface> > get_boundaries();
        template<typename Surface>
        std::shared_ptr<Surface> get_boundary(int tag); 

        std::set<std::pair<int,int>> get_patches();
        std::set<int>  get_subdomains(); // TODO: rename to subdomain tags, and write number of subdomains with retrun int 
        std::set<int> get_curves();

        Polylines& get_features() {return features; }  
        Polylines& get_borders() {return borders; }

        int number_of_subdomains(){return get_subdomains().size() ;}
        int number_of_curves(){return get_curves().size() ;}
        int number_of_patches(){return get_patches().size() ;} 
        int number_of_cells(){return c3t3.number_of_cells();}
        int number_of_surfaces(){return v.size();}

        void remove_subdomain(std::vector<int> tags);
        void remove_subdomain(int tag);              // tags 

        void add_border(Polyline_3 polyline) { borders.push_back(polyline);} 
        void add_feature(Polyline_3 polyline){ features.push_back(polyline);} 

        
        void lloyd(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true);
        void odt(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true);  
        void exude(double time_limit = 0, double sliver_bound = 0 ){ CGAL::exude_mesh_3(c3t3, sliver_bound=sliver_bound, time_limit=time_limit);} 
        void perturb(double time_limit=0, double sliver_bound=0){CGAL::perturb_mesh_3 ( c3t3, *domain_ptr.get(), time_limit=time_limit, sliver_bound=sliver_bound) ;} 


    private :
        Function_vector v; 
        std::unique_ptr<DefaultMap> map_ptr;
        std::unique_ptr<Mesh_domain> domain_ptr;
        Minimum_sphere<Kernel> min_sphere; 
        C3t3 c3t3;
        Polylines borders; 
        Polylines features;

};

void Domain::set_borders()
{
  if (this->borders.size()>0)
  { 
     domain_ptr.get()->add_features(this->borders.begin(), this->borders.end());
  }
}
void Domain::set_features()
{
  if (this->borders.size()>0)
  { 
     domain_ptr.get()->add_features(this->features.begin(), this->features.end());
  }
}

std::set<int> Domain::get_curves()
{
    const Tr& tr = c3t3.triangulation();
    std::set<int> result;
    for( auto eit = tr.finite_edges_begin(); eit != tr.finite_edges_end(); ++eit) 
    {   
         result.insert(static_cast<int>(c3t3.curve_index(*eit))); 
    } 
    return result;       

}

template<typename Surface>
Domain::Domain( std::vector<Surface> surfaces )
{
    for(typename std::vector<Surface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       Polyhedron polyhedron;
       sit->get_polyhedron(polyhedron);
       min_sphere.add_polyhedron(polyhedron);
       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       this->v.push_back(polyhedral_domain);
    }
    map_ptr = std::unique_ptr<DefaultMap>( new  DefaultMap()) ; 

    Function_wrapper wrapper(this->v, *map_ptr.get());
    domain_ptr = std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));

}
template<typename Surface>
Domain::Domain( std::vector<Surface> surfaces , AbstractMap& map )
{
 
    for(typename std::vector<Surface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       Polyhedron polyhedron; 
       sit->get_polyhedron(polyhedron);
       min_sphere.add_polyhedron(polyhedron);
       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       this->v.push_back(polyhedral_domain);
    }

    Function_wrapper wrapper(this->v,map);
    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));

}
template<typename Surface>
Domain::Domain(Surface &surface) 
{
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);

    min_sphere.add_polyhedron(polyhedron);
    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);


    this->v.push_back(polyhedral_domain);
    map_ptr = std::unique_ptr<DefaultMap>( new  DefaultMap());

    Function_wrapper wrapper(this->v, *map_ptr.get());

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox())); 


}

void Domain::create_mesh(double cell_size, double facet_size,double facet_angle,  double facet_distance,double cell_radius_edge_ratio)
{

    set_borders();
    set_features();

    std::cout << "Start meshing" << std::endl;




    Mesh_criteria criteria(CGAL::parameters::facet_angle=facet_angle ,
                           CGAL::parameters::facet_size =facet_size,
                           CGAL::parameters::facet_distance=facet_distance,
                           CGAL::parameters::cell_radius_edge_ratio=cell_radius_edge_ratio,
                           CGAL::parameters::cell_size=cell_size );


    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria,CGAL::parameters::no_exude()); 


    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();
    std::cout << "Done meshing" << std::endl;

}



void Domain::create_mesh(const double mesh_resolution )
{
    
    set_borders();
    set_features();

    double r = min_sphere.get_bounding_sphere_radius(); 
    const double cell_size = r/mesh_resolution;
    std::cout << "Cell size: " << cell_size << std::endl;

    Mesh_criteria criteria(CGAL::parameters::edge_size = cell_size,
                                       CGAL::parameters::facet_angle = 30.0,
                                       CGAL::parameters::facet_size = cell_size,
                                       CGAL::parameters::facet_distance = cell_size/10.0, 
                                       CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                       CGAL::parameters::cell_size = cell_size);

    std::cout << "Start meshing" << std::endl;
    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria,CGAL::parameters::no_exude());
   

    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();
    std::cout << "Done meshing" << std::endl;
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
/*
template<typename Surface>
void Domain::add_surface_points(Surface &surface) 
{ 

   typedef Tr::Bare_point Bare_point;
   typedef std::vector<std::pair<Bare_point, int> > Initial_points_vector;
   typedef typename Surface::Index SM_Index;
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
}*/


void Domain::remove_subdomain(int tag) 
{ 
   std::vector<int> temp;
   temp.push_back(tag);
   remove_subdomain(temp);
}




void Domain::remove_subdomain(std::vector<int> tags)
{

  int before = c3t3.number_of_cells();
  assert(before!=0);
 
  std::vector<std::tuple<Cell_handle,int,int>> rebind;

  for( auto j = tags.begin(); j!=tags.end(); ++j)
  {
    Subdomain_index temp(*j) ; 
    for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
    {
      for (std::size_t i = 0; i < 4; i++)
      { 
         if(std::find(tags.begin(), tags.end(), static_cast<int>(c3t3.subdomain_index(cit))  ) == tags.end() )          
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
  // rebind leftovers
  for ( auto cit = rebind.begin(); cit!=rebind.end(); ++cit)
  {
      Cell_handle cn = std::get<0>(*cit);
      int s =std::get<1>(*cit);
      Subdomain_index ck = Subdomain_index(std::get<2>(*cit));            
      c3t3.remove_from_complex(cn,s);
      c3t3.add_to_complex(cn,s,Surface_patch_index(ck,0) );  
  }

  int after = c3t3.number_of_cells();

  std::cout << "Number of removed subdomain cells : " << (before -after) << std::endl;
  c3t3.rescan_after_load_of_triangulation(); 
}



void Domain::add_sharp_border_edges(Polyhedron& polyhedron, double threshold)
{ 
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
template<typename Surface>
void Domain::add_sharp_border_edges(Surface& surface, double threshold) 
{ 
  Polyhedron polyhedron;
  surface.get_polyhedron(polyhedron);
  add_sharp_border_edges(polyhedron, threshold);
}

std::set<int>  Domain::get_subdomains()
{
   std::set<int> sd_indices;
   for(Cell_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
   {
        sd_indices.insert(static_cast<int>(c3t3.subdomain_index(cit)) );
   }
   return sd_indices;
}
std::set<std::pair<int,int>>  Domain::get_patches()
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
template<typename Surface>
std::shared_ptr<Surface> Domain::get_boundary(int tag)
{
   std::vector<Face> faces;
   Polyline_3 points;
   Subdomain_index useless = Subdomain_index(tag);
   facets_in_complex_3_to_triangle_soup(c3t3, useless, points, faces, true, false);  
   std::shared_ptr<Surface> surf(new  Surface(points, faces)) ;
   return surf;
}
template<typename Surface>
std::vector<std::shared_ptr<Surface>> Domain::get_boundaries() 
{ 
   std::vector<Point_3> points;
   std::vector<Face> faces;
   std::vector<std::shared_ptr<Surface>> patches;
   for ( auto i : get_patches() ) 
   {
      points.clear(); 
      faces.clear();
      if (i.first != i.second)
      { 
         facets_in_complex_3_to_triangle_soup_(c3t3,Surface_patch_index(i.first,i.second),points,faces);
         std::shared_ptr<Surface> surf(new Surface(points,faces)); 
         patches.push_back(surf);
      }
   }

   return patches;
}

void Domain::lloyd(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze )
{CGAL::lloyd_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 

void Domain::odt(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze) 
{CGAL::odt_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 


#endif
