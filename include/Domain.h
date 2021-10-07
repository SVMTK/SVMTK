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

/* --- Includes -- */
#include "Polyhedral_vector_to_labeled_function_wrapper.h"

/* -- CGAL Bounding Volumes -- */
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>

/* -- CGAL 3D Mesh Generation-- */ 
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

/* -- CGAL Mesh_3 -- */ 
#include <CGAL/Mesh_3/Detect_polylines_in_polyhedra.h>
#include <CGAL/Mesh_3/polylines_to_protect.h>
#include <CGAL/IO/read_ply_points.h>

/**
 * @brief Transform facets with a specific tag to points and facet connections.
 *
 * 
 * Based on CGAL version with similar name, but uses facet tags instead of cell tags. 
 * @see [facets_in_complex_3_to_triangle_soup](https://github.com/CGAL/cgal/blob/master/Mesh_3/include/CGAL/IO/facets_in_complex_3_to_triangle_mesh.h)
 * 
 * @param[in] c3t3 the mesh srtucture stored in the Domain class Obejct
 * @param[in] sf_index mesh facet tag. 
 * @param[out] points vector of points 
 * @param[out] faces vector of faces, i.e. std::vector<std::size_t> with 3 elements.
 * @param[in] normals_point_outside_of_the_subdomain determines the orientation of the output faces  
 * 
 * @relatesalso SVMTK Domain class.
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

/** 
 * TODO : Improve compilation
 *  @brief Writes the stored mesh to a medit file.
 *   
 *  Based on CGAL output_to_medit, but writes more information to file.
 *  @see [output_to_medit] (https://github.com/CGAL/cgal/blob/master/Mesh_3/include/CGAL/IO/File_medit.h)  
 *  The additional information is internal facets, internal edges
 *  and edge tag. This is done so that the conversion to FEniCS mesh format 
 *  is easier.
 * 
 *  @param c3t3 the mesh structure stored in the Domain class Obejct
 *  @param vertex_pmap 
 *  @param facet_pmap
 *  @param cell_pmap
 *  @param facet_twice_pmap 
 *  @param print_each_facet_twice
 *  @param save_edges 
 *  @return none, write mesh to textfile with extension .mesh.
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
 * \struct
 *
 * @brief Used to store surface points and compute 
 *        the minimum bounding radius required to 
 *        enclose all of the added surface points. 
 */
template<typename Kernel>
struct Minimum_sphere
{
           typedef typename CGAL::Min_sphere_of_spheres_d_traits_3<Kernel, typename Kernel::FT> MinSphereTraits;
           typedef typename CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
           typedef typename MinSphereTraits::Sphere Sphere;

           /**
            * @brief Adds surface points to the struct  
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
            * @brief Computes the minimum bounding radius required to enclose the added surface points
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
 * \class Domain
 * The SVMTK Domain class is used to create and tetrahedra mesh in 3D. 
 * For auxilary purpose, most of the CGAL declarations are done inside the Domain class.
 *
 *
 *
 *
 *
 * The class utilize the CGAL Exact predicates inexact constructions kernel.
 * @see(CGAL::Exact_predicates_inexact_constructions_kernel)[https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html] 
 * 
 * SVMTK Domain class uses a nested set of CGAL mesh domain in order to get
 * the required properties of the tetrahedra mesh. These properties include 
 * user specific subdomain labels, and the properties of 1D features, such as 
 * sharp edges.
 *
 * @see  (Mesh_domain_with_polyline_features_3)[https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Mesh__domain__with__polyline__features__3.html]
 *       (Polyhedral_mesh_domain_with_features_3)[https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Polyhedral__mesh__domain__with__features__3.html]
 *       (Labeled_Mesh_Domain)[https://doc.cgal.org/latest/Mesh_3/classCGAL_1_1Labeled__mesh__domain__3.html]
 *   
 * 
 * Mesh_Criteria 
 * @see (CGAL:: [https://doc.cgal.org/latest/Mesh_3/classMeshCriteria__3.html]
 *
 * Triganulation data structure 
 *@see ( CGAL::Mesh_triangulation_3)[https://doc.cgal.org/latest/Triangulation_3/index.html]
 *
 * Mesh_complex_3_in_triangulation_3
 *@see (CGAL::Mesh_complex_3_in_triangulation_3)[https://doc.cgal.org/latest/Mesh_3/classMeshComplex__3InTriangulation__3.html]
 *
 *    
 */
class Domain {
    public :
        // Overview of CGAL classe definitions and auxiliary definitions.
        // See (https://doc.cgal.org/latest/Kernel_23/group__kernel__predef.html)
        
        typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

        typedef Kernel::Point_3 Point_3;
        typedef Kernel::FT FT;
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
        // Experimental
        typedef CGAL::Mesh_constant_domain_field_3<Tr::Geom_traits,
                                          Mesh_domain::Index> Sizing_field;
        //
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
        Domain(Surface& surface,double error=1.e-7);
        template<typename Surface>
        Domain(std::vector<Surface> surfaces,double error=1.e-7);
        template<typename Surface>
        Domain(std::vector<Surface> surfaces, std::shared_ptr<AbstractMap> map,  double error=1.e-7);

        ~Domain() { for( auto vit : this->v){delete vit;}v.clear();}        

        void create_mesh(const double mesh_resolution );
        void create_mesh(double edge_size,double cell_size, double facet_size,double facet_angle,  double facet_distance,double cell_radius_edge_ratio);
        void save(std::string outpath, bool save_1Dfeatures); 

        /**
        * @breif Returns the minimum bounding sphere for all added surfaces in the constructor 
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
        std::shared_ptr<Surface> get_interface(std::pair<int,int> interface );

        template<typename Surface>
        std::shared_ptr<Surface> get_boundary(int tag); 

        std::vector<std::pair<int,int>> get_patches();
        std::set<int>  get_subdomains(); 
        std::set<int> get_curves();

        // FIXME
       template<typename Surface, typename Plane_3>
       void add_sharp_border_edges(Surface& surface,Plane_3 plane); // add const
       
        /**
         * 
         * The following functions returns an integer of a query.
         * number_of_subdomains() returns number of different cell tags 
         * number_of_curves() returns number of different edge tags 
         * number_of_patches() returns number of different facet tags 
         * number_of_surfaces() returns number surfaces used as input    
         * number_of_cells() returns number of cells in 3D mesh object c3t3;
         * number_of_facets() returns number of facets in 3D mesh object c3t3;
         * number_of_vertices() returns number of vertices in the triangulation of 
                                the 3D mesh object c3t3;        
         */
        int number_of_subdomains(){return get_subdomains().size() ;}
        int number_of_curves(){return get_curves().size() ;}
        int number_of_patches(){return get_patches().size() ;} 
        int number_of_cells(){return c3t3.number_of_cells();}
        int number_of_facets(){return c3t3.number_of_facets();}
        int number_of_vertices(){return c3t3.triangulation().number_of_vertices();}       
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
         * @brief CGAL function for excude optimazation of the constructed mesh.  
         * @see  (excude_optimize_mesh)[https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html]
         * @param time_limit 
         * @param sliver_bound 
         * @return none 
         */
        void exude(double time_limit = 0, double sliver_bound = 0 ){ assert_non_empty_mesh_object(); CGAL::exude_mesh_3(c3t3, sliver_bound=sliver_bound, time_limit=time_limit);} 
        
        /**
         * @brief CGAL function for perturb optimazation of the constructed mesh.  
         * @see  (perturb_mesh)[https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html]
         * @param time_limit 
         * @param sliver_bound 
         * @return none 
         */
        void perturb(double time_limit=0, double sliver_bound=0){assert_non_empty_mesh_object(); CGAL::perturb_mesh_3 ( c3t3, *domain_ptr.get(), time_limit=time_limit, sliver_bound=sliver_bound) ;} 

        void protect_borders();

        bool assert_non_empty_mesh_object();
       
    
        void rebind_missing_facets();
 
        template <typename Surface>
        void boundary_segmentations(int subdomain_tag,double angle_in_degree);
        template <typename Surface>
        void boundary_segmentations(double angle_in_degree);
        
        template <typename Surface>
        void boundary_segmentations(std::pair<int,int> interface, double angle_in_degree);
        
        
        std::vector<double> dihedral_angles();
        std::vector<double> radius_ratio();
        
        
        /**
         * @brief Returns the maximum and minimum dihedral angles 
         * of cells in mesh.
         *  
         * @param none. 
         * @return a pair of double with first as the minimum and second 
         * as maximum 
         */  
        std::pair<double,double > dihedral_angles_min_max()
        { 
          auto values = dihedral_angles();     
          
          double min = *std::min_element(values.begin(),values.end());
          double max = *std::max_element(values.begin(),values.end());    
          return std::make_pair(min,max);
        }
                
        /**
         * @brief Returns the maximum and minimum radius ratio 
         * of cells in mesh.
         *  
         * @param none. 
         * @return a pair of double with first as the minimum and second 
         * as maximum 
         */   
        std::pair<double,double > radius_ratio_min_max()
        { 
          auto values = radius_ratio();          
          double min = *std::min_element(values.begin(),values.end());
          double max = *std::max_element(values.begin(),values.end());    
          return std::make_pair(min,max);
          
        }
        int remove_isolated_vertices(bool remove_domain=false);

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
 * @brief Removes vertices that is not connected to any cells in the mesh.
 *
 * A feature in CGAL is that not all vetices in triangulation is part 
 * of the mesh complex. 
 *
 * @param c3t3 the mesh structure stored in the Domain class Obejct
 * @param remove_domain indication if a subdomain is removed, avoiding warning. 
 * @return integer number of vertices removed.
 */

inline int Domain::remove_isolated_vertices(bool remove_domain)
{ 

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
  double vertices_removed_ratio = 1.0 - (double)after/(double)before ;
  
  if (vertices_removed_ratio > 0.01 and !remove_domain)   
  {
       std::cout<<"There were a significant number of isolated vertices, and the user should inspect the mesh."<< std::endl;  
       std::cout<<"Methods to decrease number of isolated vertices : \n \t 1.isotropic remeshing \n \t 2.increase the mesh resolution or \n \t 3. specific mesh parameters to decrease number of isolated vertices."<< std::endl;  
  }
  return (before - after);
}


/**
 * @brief Computes the dihedral angle for all tetrahedron cells in complex (c3t3).
 * @param none 
 * @return vector of dihedral angles in degree
 */
inline std::vector<double> Domain::dihedral_angles()
{
  assert_non_empty_mesh_object();
  const Tr& tr = c3t3.triangulation();
  std::vector<double> result;
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
      result.push_back(static_cast<double>(CGAL::Mesh_3::minimum_dihedral_angle(tr.tetrahedron(cit),Kernel())));

  }
  return result;
}

/**
 * @brief Computes the radius ratio for all tetrahedron cells in complex (c3t3).
 * @param none 
 * @return vector of radius ratios.
 */
inline std::vector<double> Domain::radius_ratio()
{
  assert_non_empty_mesh_object();
  const Tr& tr = c3t3.triangulation();
  std::vector<double> result;
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
     result.push_back(static_cast<double>(CGAL::Mesh_3::radius_ratio(tr.tetrahedron(cit),Kernel() )));
  }
  return result;
}


/**
 * @brief Rebinds missing facets. 
 * 
 * The mesh may occasionally not have all the facets added in the complex.
 * This function is used to add these facets to the complex. It is important to note 
 * that this bug will only affect functions using facets_in_complex, i.e. no impact
 * regarding writing to file
 * 
 * @param none template<typename Surface>
 * @bool return true if none empty, and false if empty.
 */
inline void Domain::rebind_missing_facets()
{
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
}

/**
 * @brief Checks if 3D mesh object c3t3 is empty.
 * If empty, the function will throw an error.
 * @param none 
 * @bool return true if none empty, and false if empty.
 *
 * @throws EmptyMeshError if mesh is empty.
 */
inline bool Domain::assert_non_empty_mesh_object()
{
     if ((number_of_cells()+number_of_facets()+number_of_vertices())==0)
     {  throw  EmptyMeshError("3D mesh object is empty.");
        return false;
     }
     else
        return true;
}


/**
 * @brief Segments the boundary of a specified subdomain tag.
 * 
 * Calls Surface::surface_segmentation on the boundary surfaces specifiec by input, and 
 * updates the boundary facets with the segementation tags. 
 * 
 * @tparam SVMTK Surface object 
 * @param subdmain_tag used to obtain the boundary of subdomain with tag.
 * @param angle_in_degree the threshold angle used to detect sharp edges.
 * @return none, updates surface_patches of mesh.  
 * FIXME
 */
 
template <typename Surface>
inline void Domain::boundary_segmentations(std::pair<int,int> interface, double angle_in_degree)
{

  assert_non_empty_mesh_object();
  
  const Tr& tr = c3t3.triangulation();
  
  auto patches = get_patches();

  auto p = std::minmax_element(patches.begin(),patches.end()); //subdomain+1?

  int tag_ = 1+p.second->first;

  std::shared_ptr<Surface> surf;
  
  
  surf =  get_interface<Surface>(interface); 
  
  
    
    
  surf.get()->set_outward_face_orientation();
  surf.get()->save("jfam.off");
  surf.get()->fill_holes();
  surf.get()->save("jfam2.off");
  auto T2tagmap =surf.get()->surface_segmentation(tag_,angle_in_degree); // error

  int i,j,k,n;
  
  Cell_handle ch,cn; 
  Vertex_handle vh1,vh2,vh3;
  
  std::map<int,int> temp;
  
  for ( auto pit : T2tagmap) 
  {  
     Weighted_point wp1(pit.first[1]); 
     Weighted_point wp2(pit.first[2]); 
     Weighted_point wp3(pit.first[3]); 
 
     tr.is_vertex(wp1,vh1);
     tr.is_vertex(wp2,vh2);
     tr.is_vertex(wp3,vh3);    
               
     if ( tr.is_facet(vh1,vh2,vh3,ch,i,j,k) ) 
     {      
     n = (6- (i+j+k) ); // 3+2+1 -> 0 , 3+1+0 -> 2 etc.
     
     cn = ch->neighbor(n);
     
     //Subdomain_index ci = static_cast<int>(c3t3.subdomain_index(ch));     
     //Subdomain_index cj = static_cast<int>(c3t3.subdomain_index(cn)); 
     
     //if (!c3t3.is_in_complex(ch,n))
     //   continue;
     

     //if ( ( ci!=interface.first and cj!=interface.second ) or (cj!=interface.first and ci!=interface.second ) )
     //   continue;
     
     

     c3t3.remove_from_complex(ch,n);
     c3t3.add_to_complex(ch, n, Surface_patch_index(pit.second.first,pit.second.second) ); 
     }
     /*if( ci!=cj )
     {   
        if ( cj==0 ) 
        {
          
           c3t3.remove_from_complex(ch,n); 
           c3t3.add_to_complex(ch, n, Surface_patch_index(pit.second.first,pit.second.second)  ); 
      
        }   
        else if ( ci==0)
        {
           
           c3t3.remove_from_complex(ch,n);template<typename Surface>
           c3t3.add_to_complex(ch, n, Surface_patch_index(pit.second.first,pit.second.second) ); 

        }
        
    }*/
     
 
  }


  return;
}



/**
 * @brief Segments the boundary of a specified subdomain tag.
 * 
 * Calls Surface::surface_segmentation on the boundary surfaces specifiec by input, and 
 * updates the boundary facets with the segementation tags. 
 * 
 * @tparam SVMTK Surface object 
 * @param subdmain_tag used to obtain the boundary of subdomain with tag.
 * @param angle_in_degree the threshold angle used to detect sharp edges.
 * @return none, updates surface_patches of mesh.  
 * 
 */
template <typename Surface>
inline void Domain::boundary_segmentations(int subdomain_tag, double angle_in_degree)
{

  assert_non_empty_mesh_object();
  
  const Tr& tr = c3t3.triangulation();
  
  auto patches = get_patches();


  auto p = std::minmax_element(patches.begin(),patches.end()); 

  int tag_ = 1+p.second->first;
  
  std::shared_ptr<Surface> surf;
  
  surf =  get_boundary<Surface>(subdomain_tag);
 
  
  auto T2tagmap =surf.get()->surface_segmentation(tag_,angle_in_degree); // error

  int i,j,k,n;
  
  Cell_handle ch,cn; 
  Vertex_handle vh1,vh2,vh3;
  
  for ( auto pit : T2tagmap) 
  {  
     Weighted_point wp1(pit.first[1]); 
     Weighted_point wp2(pit.first[2]); 
     Weighted_point wp3(pit.first[3]); 
 
     tr.is_vertex(wp1,vh1);
     tr.is_vertex(wp2,vh2);
     tr.is_vertex(wp3,vh3);    
               
     tr.is_facet(vh1,vh2,vh3,ch,i,j,k);
          
     n = (6- (i+j+k) ); // 3+2+1 -> 0 , 3+1+0 -> 2 etc.
     
     cn = ch->neighbor(n);
     Subdomain_index ci = static_cast<int>(c3t3.subdomain_index(ch));     
     Subdomain_index cj = static_cast<int>(c3t3.subdomain_index(cn)); 


     if( ci!=cj )
     {   
        if ( cj==0 ) 
        {           
           c3t3.remove_from_complex(ch,n); 
           c3t3.add_to_complex(ch, n, Surface_patch_index(pit.second.first,pit.second.second)  ); 
           
        }   
        else if ( ci==0)
        {
           c3t3.remove_from_complex(ch,n);
           c3t3.add_to_complex(ch, n, Surface_patch_index(pit.second.first,pit.second.second) ); 

        }
        
    } 
     
 
  }


  return;
}


/**
 * @brief Segments the boundary of the stored mesh.
 * 
 * Calls Surface::surface_segmentation on the boundary surfaces specifiec by input, and 
 * updates the boundary facets  with the segementation tags 
 * 
 * @tparam SVMTK Surface object 
 * @param angle_in_degree the threshold angle used to detect sharp edges.
 * @return none, updates surface_patches of mesh.  
 * 
 */
template <typename Surface>
inline void Domain::boundary_segmentations(double angle_in_degree)
{

  auto tags = get_subdomains();
  for ( auto tag : tags ) 
  {
     boundary_segmentations<Surface>(tag,angle_in_degree);
  }
  
}


/**
 * @brief Sets 1-D features located on the boundary. This is done automatically 
 * in the function create_mesh.
 *
 * @param none 
 * @return void 
 */
inline void Domain::set_borders()
{
  if (this->borders.size()>0)
  { 
     protect_borders();
     domain_ptr.get()->add_features(this->borders.begin(), this->borders.end());
  }
}

/**
 * @brief Sets 1-D features located not on the boundary. This is done automatically 
 * in the function create_mesh.
 * @param none template<typename Surface>
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
 * TODO : rename 
 * @brief Returns a set of integer of the curve tags in the triangulation.
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
 * 
 * @param surface a vector of SVMTK class Surface objects defined in Surface.h 
 */
template<typename Surface>
Domain::Domain( std::vector<Surface> surfaces ,double error)
{
    for(typename std::vector<Surface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       if (!sit->does_bound_a_volume())
          sit->fill_holes();
       
       
       Polyhedron polyhedron;
       sit->get_polyhedron(polyhedron);
       min_sphere.add_polyhedron(polyhedron);
       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       this->v.push_back(polyhedral_domain);
    }
    map_ptr = std::shared_ptr<DefaultMap>( new  DefaultMap()) ; 

    Function_wrapper wrapper(this->v, map_ptr);

    domain_ptr=std::unique_ptr<Mesh_domain>( new Mesh_domain(
               Labeled_Mesh_Domain(wrapper,wrapper.bbox(),FT(error) ))); 
}

/**
 *
 * @param surfaces a vector of SVMTK class Surface objects defined in local header Surface.h 
 * @param map a smart pointer to SVMTK virtuell class AbstractMap defind in local header SubdomainMap.h
 */
template<typename Surface>
Domain::Domain( std::vector<Surface> surfaces , std::shared_ptr<AbstractMap> map, double error )
{
 
    for(typename std::vector<Surface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       if (!sit->does_bound_a_volume())
          sit->fill_holes();
        
       
       Polyhedron polyhedron; 
       sit->get_polyhedron(polyhedron);
       min_sphere.add_polyhedron(polyhedron);
       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       this->v.push_back(polyhedral_domain);
    }

    map_ptr = std::move(map);
    Function_wrapper wrapper(this->v,map_ptr);

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain( Labeled_Mesh_Domain(wrapper,wrapper.bbox(),FT(error) ))); 
}

/**
 *
 * @param surface SVMTK class Surface object defined in local header Surface.h
 */
template<typename Surface>
Domain::Domain(Surface &surface,double error) 
{
    if (!surface.does_bound_a_volume())
        surface.fill_holes();
       
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);

    min_sphere.add_polyhedron(polyhedron);
    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);

    this->v.push_back(polyhedral_domain);
    map_ptr = std::shared_ptr<DefaultMap>( new  DefaultMap());

    Function_wrapper wrapper(this->v,map_ptr);
    //

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain( Labeled_Mesh_Domain(wrapper,wrapper.bbox(),FT(error) ))); 
}

/** 
 * @brief Creates the mesh stored in the class member variable c3t3. 
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
    
    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria,CGAL::parameters::no_exude(),  
                                                               CGAL::parameters::no_perturb(),
                                                               CGAL::parameters::features(),
                                                               CGAL::parameters::non_manifold()); 
    
    remove_isolated_vertices();
    c3t3.rescan_after_load_of_triangulation();
    rebind_missing_facets();        
    std::cout << "Done meshing" << std::endl;
}

/**
 *  @brief Creates the mesh stored in the class member variable c3t3. 
 *  @note The mesh criteria is set based on mesh_resolution and the minimum bounding radius of the mesh.
 *
 *  @param mesh_resolution a value determined 
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
   
    remove_isolated_vertices();
    c3t3.rescan_after_load_of_triangulation();
    rebind_missing_facets();
    std::cout << "Done meshing" << std::endl;

}

/**
 * @brief Writes the mesh stored in the class member variable c3t3 to file.
 * 
 * The interface tags are loaded from SubdomainMap added in the constructor. 
 * If there are no interfaces in SubDomainMap, then default interfaces are 
 * selected.
 *
 * @param outpath the path to the output file.
 * @param save_1Dfeatures option to save the edges with tags.
 */
inline void Domain::save(std::string outpath,bool save_1Dfeatures)
{
    assert_non_empty_mesh_object();
    
    std::ofstream  medit_file(outpath);
    typedef CGAL::Mesh_3::Medit_pmap_generator<C3t3,false,false> Generator;
    typedef Generator::Cell_pmap Cell_pmap;
    typedef Generator::Facet_pmap Facet_pmap;
    typedef Generator::Facet_pmap_twice Facet_pmap_twice;
    typedef Generator::Vertex_pmap Vertex_pmap;
 
    Cell_pmap cell_pmap(c3t3);
    Facet_pmap facet_pmap(c3t3, cell_pmap); 
    Facet_pmap_twice  facet_twice_pmap(c3t3,cell_pmap);
    Vertex_pmap vertex_pmap(c3t3, cell_pmap,facet_pmap); // 

    std::map<std::pair<int,int>,int> facet_map = this->map_ptr->make_interfaces(this->get_patches());
    

    
    output_to_medit_(medit_file,c3t3, vertex_pmap, facet_map, cell_pmap, facet_twice_pmap , false, save_1Dfeatures) ;
    medit_file.close();
}

/**
 * @brief Removes all cells in the mesh with a specified integer tag , but perserves the 
 * interface tags as if no cells were removed.  
 * 
 * @param tag removes cells with this integer tag
 * @return void.
 * @overload
 */
inline void Domain::remove_subdomain(int tag) 
{ 
   assert_non_empty_mesh_object();

   std::vector<int> temp;
   temp.push_back(tag);
   remove_subdomain(temp);
}

/**
 * @brief Removes all cells in the mesh with tags in a vector, but perserves the 
 * interface tags as if no cells were removed.
 * @param tags vector of cell tag to be removed. 
 * @return void
 * @overload
 */
inline void Domain::remove_subdomain(std::vector<int> tags)
{
  assert_non_empty_mesh_object();
  
  int before = c3t3.number_of_cells();
  std::vector<std::tuple<Cell_handle,int,int,int>> rebind;  
  int spindex1,spindex2;
  // Iterate over subdomains to be removed, and stores connected facets (Cell_handle,int) and patches (int,int)
  for( auto j = tags.begin(); j!=tags.end(); ++j)
  {
    Subdomain_index temp(*j) ; 
    // iterate over all cells, intended to be all other cells.
    for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
    {
      for (std::size_t i = 0; i < 4; i++)
      { 
         spindex1 = static_cast<int>(c3t3.surface_patch_index(cit,i).first);
         spindex2 = static_cast<int>(c3t3.surface_patch_index(cit,i).second); 
        
         // True if the cell is not a cell to be removed 
         if(std::find(tags.begin(), tags.end(), static_cast<int>(c3t3.subdomain_index(cit))  ) == tags.end() and c3t3.is_in_complex(cit))           
         {
           // True if neighbor cell is a cell to be removed
           if(c3t3.subdomain_index(cit->neighbor(i))==temp or !c3t3.is_in_complex(cit->neighbor(i)))
           {     
               rebind.push_back(std::make_tuple(cit, i,spindex1,spindex2));               
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
  remove_isolated_vertices(true);

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
 *
 * @brief adds sharp border edges of polyhedron to mesh.
 *
 * Checks polyhedron for sharp edges, if these edges do not conflict/intersect
 * previously stored edges, then the edges are stored in this->borders.
 *   
 * @note The use of 1D features in combination with ill-posed meshing paramteres can cause segmentation fault (crash). 
 *
 * @param polyhedron triangulated surface mesh in 3D (@see Domain::polyhedron).
 * @param threshold that determines sharp edges.
 * @precondition edges can not intersect in a non conforming way.
 * @overload
 */
inline void Domain::add_sharp_border_edges(Polyhedron& polyhedron, double threshold)
{ 
 typedef boost::property_map<Polyhedron, CGAL::edge_is_feature_t>::type EIF_map;

 EIF_map eif = get(CGAL::edge_is_feature, polyhedron);
 
 CGAL::Polygon_mesh_processing::detect_sharp_edges(polyhedron,threshold, eif); 

 Point_3  p1,p2;
 Polylines temp;

 bool intersecting=false;
 for(boost::graph_traits<Polyhedron>::edge_descriptor e : edges(polyhedron))
 {
    if(get(eif, e))
    {
       Polyline_3 polyline;

       p1 = source(e,polyhedron)->point();
       p2 = target(e,polyhedron)->point();
      
       polyline={p1,p2};
       intersecting=false;
       for ( auto pline : this->borders ) 
       {
           if (CGAL::Polygon_mesh_processing::do_intersect(polyline,pline)) 
              intersecting=true;      
       } 
       if (!intersecting)
          temp.push_back(polyline);
    
    }
 }

 if (temp.size()==0)
   std::cout <<"Warning, new edges intersects with existing edges."<<std::endl;
 else
   this->borders.insert(this->borders.end(),temp.begin(),temp.end());    
}




/**
 * @brief Adds sharp border edges from a SVMTK Surface object.
 *
 * @tparam SVMTK Surface class
 * @param surface Surface object defined in Surface.h .
 * @param threshold angle in degree (0-90) bewteen to connected edges.
 * @overload
 */
template<typename Surface>
void Domain::add_sharp_border_edges(Surface& surface, double threshold) // add const
{ 

  surface.collapse_edges();
  Polyhedron polyhedron;  
  surface.get_polyhedron(polyhedron);
  add_sharp_border_edges(polyhedron, threshold);
}



/**
 * @brief Adds sharp border edges from a SVMTK Surface object.
 *
 * @tparam SVMTK Surface class
 * @param surface Surface object defined in Surface.h .
 * @param threshold angle in degree (0-90) bewteen to connected edges.
 * @overload
 */
template<typename Surface, typename Plane_3>
void Domain::add_sharp_border_edges(Surface& surface,Plane_3 plane) // add const
{ 

  surface.collapse_edges(); // FIXME
  
  auto polylines = surface.polylines_in_plane(plane);
  for ( auto polyline : polylines) 
  {
    add_feature(polyline);
  } 
  
}


/**
 * @brief Combines borders edges to a connected polyline.
 *
 * If the number polylines remains the same, 
 * then the borders are not updated, this is to 
 * avoid core dump.  
 * 
 * @param none 
 * @return void updates borders variable 
 */
inline void Domain::protect_borders() 
{
   Polylines temp;
   polylines_to_protect(temp, this->borders.begin() ,this->borders.end());
   if ( temp.size() == this->borders.size())
      return;
   else 
      this->borders = temp;
}

/**
 * @brief Returns a set of integer that represents the cell tags in the mesh.   
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
 * @brief Returns a set of integer pairs that represents the surface facet tags in the mesh. 
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

         if (dummy.insert(std::pair<int,int>(static_cast<int>(spi.first) , static_cast<int>(spi.second) )).second and spi.first!=spi.second )
          {   sf_indices.push_back(std::pair<int,int>(static_cast<int>(spi.first) , static_cast<int>(spi.second) ));

          }
      }
   }
   std::sort(sf_indices.begin() ,sf_indices.end()) ;
   return sf_indices;
}

// 
/**
 * @brief Returns the boundary of a subdomain tag stored in a Surface object.
 *
 * @tparam SVMTK Surface object   
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
 * @brief Iterates over all surface boundaries of subdomains and stores and returns it as a vector of Surface objects.   
 *
 * @tparam SVMTK Surface object    
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


template<typename Surface>
std::shared_ptr<Surface> Domain::get_interface(std::pair<int,int> interface ) 
{ 
   
   std::vector<Point_3> points;
   std::vector<Face> faces;
   std::vector<std::shared_ptr<Surface>> patches;

   if (interface.first == interface.second)
       throw InvalidArgumentError("There are no interfaces between similar tags.") ;

   facets_in_complex_3_to_triangle_soup_(c3t3,Surface_patch_index(interface.first,interface.second),points,faces);
   std::shared_ptr<Surface> surf(new Surface(points,faces)); 

   return surf;
}



/**
 * @brief CGAL function for lloyd optimazation of the constructed mesh. 
 *
 * @see (lloyd_optimize_mesh)[https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html]
 * @note Do not use lloyd after excude optimazation, since it will cause failure.   
 * @param time_limit 
 * @param max_iteration_number   
 * @param convergence 
 * @param freeze_bound 
 * @param do_freeze
 * @return none
 */
inline void Domain::lloyd(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze )
{   assert_non_empty_mesh_object();

CGAL::lloyd_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 

/**
 * @brief CGAL function for odt optimazation of the constructed mesh.  
 * @see  (odt_optimize_mesh)[https://doc.cgal.org/latest/Mesh_3/group__PkgMesh3Functions.html]
 *
 * @param time_limit 
 * @param max_iteration_number   
 * @param convergence 
 * @param freeze_bound 
 * @param do_freeze
 * @return none 

 */
inline void Domain::odt(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze) 
{   assert_non_empty_mesh_object();

CGAL::odt_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 

#endif
