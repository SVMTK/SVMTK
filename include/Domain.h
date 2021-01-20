// Copyright (C) 2018-2021 Lars Magnus Valnes 
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

#include <CGAL/Mesh_3/Detect_polylines_in_polyhedra.h>


#include <CGAL/Mesh_3/polylines_to_protect.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/Triangulation_3.h>

/**
 *  Based on CGAL version with similar name, but uses facet tags instead of cell tags. 
 *  Outputs the vertices and faces of a facet tag in the mesh.
 *  @param[in] c3t3 the mesh srtucture stored in the Domain class Obejct
 *  @param[in] sf_index mesh facet tag. 
 *  @param[out] points 
 *  @param[out] faces 
 *  @param[in] normals_point_outside_of_the_subdomain determines the orientation of the output faces  
 *  @see [facets_in_complex_3_to_triangle_soup](https://github.com/CGAL/cgal/blob/master/Mesh_3/include/CGAL/IO/facets_in_complex_3_to_triangle_mesh.h)
 */
template<class C3T3, class PointContainer, class FaceContainer>
void facets_in_complex_3_to_triangle_soup_(const C3T3& c3t3,
                                          const typename C3T3::Surface_patch_index sf_index,
                                          PointContainer& points,
                                          FaceContainer& faces,
                                          const bool normals_point_outside_of_the_subdomain = false)
                                       
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

/** Based on CGAL output_to_medit, but writes more information to file.
 *  The additional information is internal facets, internal edges
 *  and edge tag. This is done so that the conversion to FEniCS mesh format 
 *  is easier.
 *  @param c3t3 the mesh structure stored in the Domain class Obejct
 *  @param vertex_pmap
 *  @param facet_pmap
 *  @param cell_pmap
 *  @param facet_twice_pmap 
 *  @param print_each_facet_twice
 *  @param save_edges 
 *  @return none, write mesh to textfile with extension .mesh.
 *  @see [output_to_medit] (https://github.com/CGAL/cgal/blob/master/Mesh_3/include/CGAL/IO/File_medit.h)  
*/
template <class C3T3,
          class Vertex_index_property_map,
          typename Facet_index_property_map,
          class Facet_index_property_map_twice,
          class Cell_index_property_map>
void output_to_medit_(std::ostream& os,
                const C3T3& c3t3,
                const Vertex_index_property_map& vertex_pmap,
                Facet_index_property_map& facet_pmap,
                const Cell_index_property_map& cell_pmap,
                const Facet_index_property_map_twice& facet_twice_pmap = Facet_index_property_map_twice(),
                const bool print_each_facet_twice = false,
                const bool save_edges = true )
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
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

  os << "Triangles\n" << number_of_triangles << '\n';
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
    
    if (c3t3.is_in_complex(fit->first) or  c3t3.is_in_complex(fit->first->neighbor(fit->second) ) )  
    {
       os << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3] << ' '; 
       Surface_patch_index spi = c3t3.surface_patch_index(*fit);
       std::pair<int,int> key(static_cast<int>(spi.first) , static_cast<int>(spi.second) );
       if ( key.second> key.first){std::swap(key.first,key.second);}
       if (facet_pmap.find(key) == facet_pmap.end() )
       {
          os << 0 << '\n'; 
       }
       else
       {
          os << facet_pmap.at(key) << '\n';
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
  os << "End\n";
}

/**
 * @struct:
 * Used to store surface points and compute 
 * the minimum bounding radius required to 
 * enclose all of the added surface points 
 */
template<typename Kernel>
struct Minimum_sphere
{
           typedef typename CGAL::Min_sphere_of_spheres_d_traits_3<Kernel, typename Kernel::FT> MinSphereTraits;
           typedef typename CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
           typedef typename MinSphereTraits::Sphere Sphere;

           /**
            * Adds surface points to the struct  
            * @param polyhedron triangulated surface structure.
            * @return none 
            */            
           template< typename MeshPolyhedron_3>  
           void add_polyhedron(const MeshPolyhedron_3 &polyhedron)
           {
                for (typename MeshPolyhedron_3::Vertex_const_iterator it=polyhedron.vertices_begin();it != polyhedron.vertices_end(); ++it)
                {
                    S.push_back(Sphere(it->point(), 0.0));
                }
            } 

           /**
            * Computes the minimum bounding radius required to enclose the added surface points
            * @param none   
            * @return the radius that encloses all added surface points  
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
 * Removes vertices that is not connected to any cells in the mesh.
 * @param c3t3 the mesh structure stored in the Domain class Obejct
 * @param remove_domain indication if a subdomain is removed, avoiding warning. 
 * @return integer number of vertices removed.
 */
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

/**
 * @Class 
 * Used to create and store the resulting tetrahedra mesh.
 */
class Domain {
    public :
        // Overview of CGAL classe definitions and auxiliary definitions.
        // See (https://doc.cgal.org/latest/Kernel_23/group__kernel__predef.html)
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
        Domain(std::vector<Surface> surfaces, std::shared_ptr<AbstractMap> map);

        ~Domain() { for( auto vit : this->v){delete vit;}v.clear();}        

        void create_mesh(const double mesh_resolution );
        void create_mesh(double edge_size,double cell_size, double facet_size,double facet_angle,  double facet_distance,double cell_radius_edge_ratio);
        void save(std::string outpath, bool save_1Dfeatures); 

        /**
        * Returns the minimum bounding sphere for all added surfaces in the constructor 
        * @param none
        * @return the minimum bounding sphere for all added surfaces in the constructor 
        */
        double get_bounding_sphere_radius(){ return min_sphere.get_bounding_sphere_radius(); }

        void add_sharp_border_edges(Polyhedron& polyhedron, double threshold);
        template<typename Surface>
        void add_sharp_border_edges(Surface& surface, double threshold=60);

        template<typename Surface>
        std::vector<std::shared_ptr<Surface> > get_boundaries();
        template<typename Surface>
        std::shared_ptr<Surface> get_boundary(int tag); 

        std::vector<std::pair<int,int>> get_patches();
        std::set<int>  get_subdomains(); 
        std::set<int> get_curves();

        /**
         * The following functions returns an integer of a query.
         * number_of_subdomains() returns number of different cell tags 
         * number_of_curves() returns number of different edge tags 
         * number_of_patches() returns number of different facet tags 
         * number_of_surfaces() returns number surfaces used as input            
         */
        int number_of_subdomains(){return get_subdomains().size() ;}
        int number_of_curves(){return get_curves().size() ;}
        int number_of_patches(){return get_patches().size() ;} 
        int number_of_cells(){return c3t3.number_of_cells();}
        int number_of_surfaces(){return v.size();}
        
        void remove_subdomain(std::vector<int> tags);
        void remove_subdomain(int tag);              

        void set_borders();
        void set_features();

        void clear_borders(){this->borders.clear();}
        void clear_features(){this->features.clear();}

        void add_border(Polyline_3 polyline) { borders.push_back(polyline);} 
        void add_feature(Polyline_3 polyline){ features.push_back(polyline);} 

        Polylines& get_features() {return features; }  
        Polylines& get_borders() {return borders; }
        
        void lloyd(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true);
        void odt(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true); 
          
        /**
         * Wrapper CGAL function for odt optimazation of the constructed mesh.  
         * @see  (excude_optimize_mesh)[https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html]
         * @param time_limit 
         * @param sliver_bound 
         * @return none 
         */
        void exude(double time_limit = 0, double sliver_bound = 0 ){ CGAL::exude_mesh_3(c3t3, sliver_bound=sliver_bound, time_limit=time_limit);} 
        /**
         * Wrapper CGAL function for odt optimazation of the constructed mesh.  
         * @see  (perturb_mesh)[https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html]
         * @param time_limit 
         * @param sliver_bound 
         * @return none 
         */
        void perturb(double time_limit=0, double sliver_bound=0){CGAL::perturb_mesh_3 ( c3t3, *domain_ptr.get(), time_limit=time_limit, sliver_bound=sliver_bound) ;} 


    private :
        Function_vector v; 
        std::shared_ptr<AbstractMap> map_ptr;
        std::unique_ptr<Mesh_domain> domain_ptr;
        Minimum_sphere<Kernel> min_sphere; 
        C3t3 c3t3;
        Polylines borders; 
        Polylines features;

};

/**
 * Sets 1-D features located on the boundary. This is done automatically 
 * in the function create_mesh.
 * @param none 
 * @return none 
 */
inline void Domain::set_borders()
{
  if (this->borders.size()>0)
  { 
     domain_ptr.get()->add_features(this->borders.begin(), this->borders.end());
  }
}

/**
 * Sets 1-D features located not on the boundary. This is done automatically 
 * in the function create_mesh.
 * @param none 
 * @return none 
 */
inline void Domain::set_features()
{
  if (this->features.size()>0)
  { 
     domain_ptr.get()->add_features(this->features.begin(), this->features.end());
  }
}

/**
 * Returns a set of integer of the curve tags in the triangulation.
 * @note This should be done after a call to create_mesh. 
 * @param none 
 * @return result a set of integers representing the curve tags.
 */
inline std::set<int> Domain::get_curves()
{
    const Tr& tr = c3t3.triangulation();
    std::set<int> result;
    for( auto eit = tr.finite_edges_begin(); eit != tr.finite_edges_end(); ++eit) 
    {   
         result.insert(static_cast<int>(c3t3.curve_index(*eit))); 
    } 
    return result;       
}

/**
 * Constructor:
 * @param surface a vector of SVMTK class Surface objects defined in Surface.h 
 */
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
    map_ptr = std::shared_ptr<DefaultMap>( new  DefaultMap()) ; 

    Function_wrapper wrapper(this->v, map_ptr);
    domain_ptr = std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));
}

/**
 * Constructor:
 *
 * @param surfaces a vector of SVMTK class Surface objects defined in local header Surface.h 
 * @param map a smart pointer to SVMTK virtuell class AbstractMap defind in local header SubdomainMap.h
 */
template<typename Surface>
Domain::Domain( std::vector<Surface> surfaces , std::shared_ptr<AbstractMap> map )
{
 
    for(typename std::vector<Surface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       Polyhedron polyhedron; 
       sit->get_polyhedron(polyhedron);
       min_sphere.add_polyhedron(polyhedron);
       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       this->v.push_back(polyhedral_domain);
    }

    map_ptr = std::move(map);
    Function_wrapper wrapper(this->v,map_ptr);
    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));
}

/**
 * Constructor:
 * 
 * @param surface SVMTK class Surface object defined in local header Surface.h
 */
template<typename Surface>
Domain::Domain(Surface &surface) 
{
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);

    min_sphere.add_polyhedron(polyhedron);
    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);

    this->v.push_back(polyhedral_domain);
    map_ptr = std::shared_ptr<DefaultMap>( new  DefaultMap());

    Function_wrapper wrapper(this->v,map_ptr);
    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox())); 
}

/**
 * 
 * Creates the mesh stored in the class member variable c3t3. 
 * Specific CGAL mesh criteria parameters.
 * @param edge_size mesh criteria for the maximum edge size 
 * @param cell_size mesh criteria for the maximum cell size  
 * @param facet_size mesh criteria for the maximum facet size 
 * @param facet_angle mesh criteria for the minimum edge size 
 * @param facet_distance mesh criteria for surface approximation  
 * @param cell_radius_edge_ratio mesh criteria for the relation between cell rddius and edge
 * @see [CGAL::Mesh_criteria](https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__criteria__3.html) 
 * @overload
 */
inline void Domain::create_mesh(double edge_size,double cell_size, double facet_size,double facet_angle,  double facet_distance,double cell_radius_edge_ratio)
{
    set_borders();
    set_features();

    Mesh_criteria criteria(CGAL::parameters::edge_size = edge_size,
                           CGAL::parameters::facet_angle=facet_angle ,
                           CGAL::parameters::facet_size =facet_size,
                           CGAL::parameters::facet_distance=facet_distance,
                           CGAL::parameters::cell_radius_edge_ratio=cell_radius_edge_ratio,
                           CGAL::parameters::cell_size=cell_size );

    std::cout << "Start meshing" << std::endl;
    
    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria,CGAL::parameters::no_exude()); 
    
    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();
    std::cout << "Done meshing" << std::endl;
}

/**
 *  Creates the mesh stored in the class member variable c3t3. 
 *  @note The mesh criteria is set based on mesh_resolution and the minimum bounding radius of the mesh.
 *  @param mesh_resolution a value determined 
 *  @see [CGAL::Mesh_criteria](https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__criteria__3.html) 
 *  @overload
 */
inline void Domain::create_mesh(const double mesh_resolution )
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

/**
 * Writes the mesh stored in the class member variable c3t3 to file. 
 * The interfaces tags are loaded from SubdomainMap added in the constructor.
 * @param outpath the path to the output file.
 * @param save_1Dfeatures option to save the edges with tags.
 */
inline void Domain::save(std::string outpath,bool save_1Dfeatures)
{
    std::ofstream  medit_file(outpath);
    typedef CGAL::Mesh_3::Medit_pmap_generator<C3t3,false,false> Generator;
    typedef typename Generator::Cell_pmap Cell_pmap;
    typedef typename Generator::Facet_pmap Facet_pmap;
    typedef typename Generator::Facet_pmap_twice Facet_pmap_twice;
    typedef typename Generator::Vertex_pmap Vertex_pmap;
 
    Cell_pmap cell_pmap(c3t3);
    Facet_pmap facet_pmap(c3t3, cell_pmap); 
    Facet_pmap_twice  facet_twice_pmap(c3t3,cell_pmap);
    Vertex_pmap vertex_pmap(c3t3, cell_pmap,facet_pmap); // 

    std::map<std::pair<int,int>,int> facet_map = this->map_ptr->get_interfaces(this->get_patches() );
    
    output_to_medit_(medit_file,c3t3, vertex_pmap, facet_map, cell_pmap, facet_twice_pmap , false, save_1Dfeatures) ;
    medit_file.close();
}

/**
 * Removes all cells in the mesh with a specified integer tag , but perserves the 
 * interface tags as if no cells were removed.  
 * @param tag removes cells with this integer tag
 * @overload
 */
inline void Domain::remove_subdomain(int tag) 
{ 
   std::vector<int> temp;
   temp.push_back(tag);
   remove_subdomain(temp);
}

/**
 * Removes all cells in the mesh with tags in a vector, but perserves the 
 * interface tags as if no cells were removed.
 * @param tags vector of cell tag to be removed. 
 * @overload
 */
inline void Domain::remove_subdomain(std::vector<int> tags)
{

  int before = c3t3.number_of_cells();

  assert(before!=0);
  std::vector<std::tuple<Cell_handle,int,int,int>> rebind;  
  // Iterate over subdomains to be removed, and stores connected facets (Cell_handle,int) and patches (int,int)
  for( auto j = tags.begin(); j!=tags.end(); ++j)
  {
    Subdomain_index temp(*j) ; 
    for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
    {
      int k = static_cast<int>(c3t3.subdomain_index(cit)) ;
      for (std::size_t i = 0; i < 4; i++)
      { 

         if(std::find(tags.begin(), tags.end(), static_cast<int>(c3t3.subdomain_index(cit))  ) == tags.end() )          
         {
           if(c3t3.subdomain_index(cit->neighbor(i))==temp)
           {         
               rebind.push_back(std::make_tuple(cit, i,*j,k)); 
           }
         }
      }
    }
  }

  // removes all cells with a subdomain tag   
  for(std::vector<int>::iterator j = tags.begin(); j!=tags.end(); ++j)
  {
      for(auto cit =  c3t3.cells_in_complex_begin(Subdomain_index(*j));cit !=  c3t3.cells_in_complex_end(); ++cit)
      {  
         c3t3.remove_from_complex(cit); 
      }     
  }

  c3t3.rescan_after_load_of_triangulation(); 
  remove_isolated_vertices(c3t3,true);

  // Go through all cells and rebind factes.
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
      Cell_handle cn = cit;
      Subdomain_index ci = c3t3.subdomain_index(cn);
      for(int i=0; i<4; ++i)
      {
        Cell_handle cm = cit->neighbor(i);
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
  // Restore the patches and facets linked to deleted cells 
  for ( auto cit = rebind.begin(); cit!=rebind.end(); ++cit)
  {
      Cell_handle cn = std::get<0>(*cit);
      int s =std::get<1>(*cit);
      Subdomain_index cf = Subdomain_index(std::get<3>(*cit));            
      Subdomain_index ck = Subdomain_index(std::get<2>(*cit)); 
      c3t3.remove_from_complex(cn,s);
      if (cf>ck) 
      {
      c3t3.add_to_complex(cn,s,Surface_patch_index(cf,ck) );  
      }
      else 
      {
      c3t3.add_to_complex(cn,s,Surface_patch_index(ck,cf) );  
      }

  }
  int after = c3t3.number_of_cells();

  std::cout << "Number of removed subdomain cells : " << (before -after) << std::endl;
  c3t3.rescan_after_load_of_triangulation(); 
}

/**
 * Checks polyhedron for sharp edges, and store this as 1-D border features.   
 * @note In combination with ill-posed meshing paramteres can cause segmentation fault (crash). 
 * @param polyhedron definition found in Domain object typedefs.
 * @param threshold that determines sharp edges.
 * @overload
 */
inline void Domain::add_sharp_border_edges(Polyhedron& polyhedron, double threshold)
{ 
 typedef boost::property_map<Polyhedron, CGAL::edge_is_feature_t>::type EIF_map;
 EIF_map eif = get(CGAL::edge_is_feature, polyhedron);
 Polylines temp;
 CGAL::Polygon_mesh_processing::detect_sharp_edges(polyhedron,threshold, eif); 
 for( Polyhedron::Edge_iterator he = polyhedron.edges_begin(); he != polyhedron.edges_end() ; ++he)
 {
   if(he->is_feature_edge() ) 
   {
       Polyline_3 polyline;
       polyline.push_back(he->vertex()->point());
       polyline.push_back(he->opposite()->vertex()->point());     
       this->borders.push_back(polyline);
   }    
 }  

}

/**
 * Checks Surface object for sharp edges, and store this as 1-D border features.   
 * @note In combination with ill-posed meshing paramteres can cause segmentation fault (crash). 
 * @param surface Surface object defined in Surface.h .
 * @param threshold angle in degree (0-90) bewteen to connected edges.
 * @overload
 */
template<typename Surface>
void Domain::add_sharp_border_edges(Surface& surface, double threshold) 
{ 
  Polyhedron polyhedron;
  surface.get_polyhedron(polyhedron);
  add_sharp_border_edges(polyhedron, threshold);
}

/**
 * Returns a set of integer that represents the cell tags in the mesh.   
 * @param none 
 * @return sd_indices a set of integers that represents the cell tags in the mesh.
 */
std::set<int>  Domain::get_subdomains()
{
   std::set<int> sd_indices;
   for(Cell_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
   {
        sd_indices.insert(static_cast<int>(c3t3.subdomain_index(cit)) );
   }
   return sd_indices;
}


/**
 * Returns a set of integer pairs that represents the surface facet tags in the mesh. 
  * @param none 
 * @return sf_indices a vector of integer pairs that represents the facet tags in the mesh.
 */
std::vector<std::pair<int,int>>  Domain::get_patches()
{

   std::vector<std::pair<int,int>> sf_indices;
   std::set<std::pair<int,int>> dummy;
   for(Cell_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
   {
      for( int i =0 ; i<4 ; ++i)
      {
         Surface_patch_index spi = c3t3.surface_patch_index(cit,i) ;
    
         if (dummy.insert(std::pair<int,int>(static_cast<int>(spi.first) , static_cast<int>(spi.second) )).second 
                 and spi.first!=spi.second )
            sf_indices.push_back(std::pair<int,int>(static_cast<int>(spi.first) , static_cast<int>(spi.second) ));

      }
   }
   std::sort(sf_indices.begin() ,sf_indices.end()) ;
   return sf_indices;
}

/**
 * Returns the boundary of a subdomain tag stored in a Surface object.   
 * @param tag the subdomain tag 
 * @return surf Surface object 
 */
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

/**
 * Iterates over all surface boundaries of subdomains and stores and returns it as a vector of Surface objects.   
 * @param none 
 * @return patches vector of Surface objects 
 */
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

/**
 * Wrapper CGAL function for lloyd optimazation of the constructed mesh. 
 * @note Do not use lloyd after excude optimazation, since it will cause failure.  
 * @see (lloyd_optimize_mesh)[https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html]
 * @param time_limit 
 * @param max_iteration_number   
 * @param convergence 
 * @param freeze_bound 
 * @param do_freeze
 * @return none
 */
inline void Domain::lloyd(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze )
{CGAL::lloyd_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 

/**
 * Wrapper CGAL function for odt optimazation of the constructed mesh.  
 * @see  (odt_optimize_mesh)[https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html]
 * @param time_limit 
 * @param max_iteration_number   
 * @param convergence 
 * @param freeze_bound 
 * @param do_freeze
 * @return none 
 * @note Do not use lloyd after excude optimazation, since it will cause failure. 
 */
inline void Domain::odt(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze) 
{CGAL::odt_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 

#endif
