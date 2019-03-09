#ifndef __CGAL_MESH_CREATOR_H


#define __CGAL_MESH_CREATOR_H



// LOCAL
#include "CGALSurface.h" 
#include "SubdomainMap.h" 
#include "Polyhedral_vector_to_labeled_function_wrapper.h"

// STD 
#include <list>
#include <fstream>
#include <memory>
 
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


#include <CGAL/Mesh_3/Dump_c3t3.h>
#include <CGAL/Mesh_3/Dump_c3t3.h>








// CURRENT PLAN :
// TODO: Consider void remove_isolated_vertices() as stand alone function or class function
// TODO: Clean up header files
// UPDATE : Lip_sizing field removed



class CGALSurface;

class CGALMeshCreator {
    public :

        //---------------------------------------------------------------
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef K::Point_3 Point_3;
        typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron; 
        typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron> Polyhedral_mesh_domain_3; // Use CGAL::Polyhedral_mesh_domain_with_features_3<K> Polyhedral_mesh_domain_3;
        typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, K  > Function_wrapper; //


        typedef Function_wrapper::Function_vector Function_vector; //
        typedef CGAL::Labeled_mesh_domain_3<K> Labeled_Mesh_Domain;
        typedef CGAL::Mesh_domain_with_polyline_features_3<Labeled_Mesh_Domain> Mesh_domain; // labeled mesh_domain function wrapper 

       
        typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

        typedef Mesh_domain::Curve_segment_index Curve_index;
        typedef Mesh_domain::Corner_index Corner_index;

        typedef CGAL::Mesh_complex_3_in_triangulation_3< Tr, Corner_index, Curve_index> C3t3;

        typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
        typedef Tr::Vertex_handle Vertex_handle;
        typedef Tr::Finite_vertices_iterator Finite_vertices_iterator;
        typedef Tr::Locate_type Locate_type;

        typedef std::vector<Point_3>  Polyline_3;
        typedef std::vector<Polyline_3> Polylines;


        typedef C3t3::Subdomain_index Subdomain_index;
        typedef C3t3::Cells_in_complex_iterator Cell_iterator;

        typedef std::map<std::string, double> Parameters;

   
        struct Minimum_sphere
        {
   
           typedef CGAL::Min_sphere_of_spheres_d_traits_3<K, K::FT> MinSphereTraits;
           typedef CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
           typedef MinSphereTraits::Sphere Sphere;

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

        CGALMeshCreator(CGALSurface& surface);
        CGALMeshCreator(std::vector<CGALSurface> surfaces);
        CGALMeshCreator(std::vector<CGALSurface> surfaces, AbstractMap& map);

        ~CGALMeshCreator() {}

        void set_parameters(std::map<std::string, double> &new_parameters); // py::dict? 
        void set_parameter(std::string key, double value);

        void create_mesh();

        void create_mesh(const double mesh_resolution );

        void default_parameters() {
            parameters["mesh_resolution"]=64.0;
            parameters["facet_angle"]    = 25.0;
            parameters["facet_size"]     = 0.1;
            parameters["facet_distance"] =  0.1;
            parameters["cell_radius_edge_ratio"] = 3.0;
            parameters["cell_size"] = 0.1;
            parameters["edge_size"] = 0.1;
        }

        void default_creating_mesh();  // Remove ?? eller default

        void save(std::string OutPath); 

        void refine_mesh(const double mesh_resolution );

        void refine_mesh();

        Polylines& get_features() {return features; }  

        void add_feature(Polyline_3 polyline) { features.push_back(polyline);} // Check 

        void set_features();



        Polylines& get_borders() {return borders; }

        void set_borders(){domain_ptr.get()->add_features(get_borders().begin(), get_borders() .end());} // -> TODO: integrate or REOMVE 

        void add_borders(Polyline_3 polyline) { borders.push_back(polyline);} 
      
        void add_sharp_border_edges(Polyhedron& polyhedron);

        void add_sharp_border_edges(CGALSurface& surface); 
 
        void reset_borders(){borders.clear();} // rename

        void remove_border(){return;}   

     
        void set_features(Polylines& polylines){ set_features( polylines.begin() , polylines.end() );} 

        template< typename InputIterator> // Wrapping 
        void set_features(InputIterator begin, InputIterator end){domain_ptr->add_features(begin, end);} 

        std::shared_ptr<CGALSurface> get_boundary();
        // Direct wrapping of CGAL
        void lloyd(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true);
        void odt(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true);
        void exude( double time_limit = 0, double sliver_bound = 0 ){ CGAL::exude_mesh_3(c3t3, sliver_bound=sliver_bound, time_limit=time_limit);} 
        void perturb( double time_limit=0, double sliver_bound=0){CGAL::perturb_mesh_3 ( c3t3, *domain_ptr.get(), time_limit=time_limit, sliver_bound=sliver_bound) ;} 


        // Specific labeling of boundary cells 
        void label_boundary_cells(int btag, int ntag); // 
        void remove_label_cells(int tag);              // tags 
        
        int number_of_cells(){return c3t3.number_of_cells();}


    private :
        std::unique_ptr<Mesh_domain> domain_ptr;
        Minimum_sphere min_sphere; 

        Polylines borders; 
        Polylines features;
        Parameters parameters;  
        C3t3 c3t3;
};



template<typename C3T3>
int remove_isolated_vertices(C3T3& c3t3)
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



  std::cout<<"Number of vertices removed: "  << before - after  << std::endl;
  if ((before - after) > 10)   
  {
       // FIXME : Better 
       std::cout<<"The number of removed vertices are substantial."<< std::endl;  
       std::cout<<"This can cause the mesh to be ill-posed."<< std::endl;  
       std::cout<<"Try isotropic remeshing or increase the mesh resolution."<< std::endl;  
  }

  return (before - after);


}

// --------------------------------
//  Constructor for CGALSurfaces
// --------------------------------
CGALMeshCreator::CGALMeshCreator( std::vector<CGALSurface> surfaces )
{
    Function_vector v;

    for(std::vector<CGALSurface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       Polyhedron *polyhedron = new Polyhedron();
       sit->get_polyhedron(*polyhedron);
       min_sphere.add_polyhedron(*polyhedron);

       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(*polyhedron);
       v.push_back(polyhedral_domain);
    }
    Function_wrapper wrapper(v);
    
    domain_ptr = std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,wrapper.bbox()));
    default_parameters();
}

CGALMeshCreator::CGALMeshCreator( std::vector<CGALSurface> surfaces , AbstractMap& map )
{
    Function_vector v;

    for(std::vector<CGALSurface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
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




CGALMeshCreator::CGALMeshCreator(CGALSurface &surface) 
{
    Function_vector v;
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);

    min_sphere.add_polyhedron(polyhedron);

    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);

    v.push_back(polyhedral_domain);
    Function_wrapper wrapper(v);



    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,wrapper.bbox())); 
    default_parameters();

}

// --------------------------------
//  
// --------------------------------

std::shared_ptr<CGALSurface> CGALMeshCreator::get_boundary()
{

   std::shared_ptr<CGALSurface> surf(new  CGALSurface()) ;
   if ( c3t3.number_of_cells()>0 )
   {
        facets_in_complex_3_to_triangle_mesh( c3t3  ,surf->get_mesh()) ;
   }
   else
   {
       std::cout <<"No mesh generated, output is an empty surface file" << std::endl; 
   }
   return surf;
}



void CGALMeshCreator::set_parameters(std::map<std::string, double> &new_parameters)
{
    for (std::map<std::string, double>::iterator pit= new_parameters.begin(); pit!=new_parameters.end(); ++pit )
    {
        parameters[pit->first] = static_cast<double>(pit->second);
    }
}
void CGALMeshCreator::set_parameter(std::string key , double value )
{
   parameters[key] = value;
}

void CGALMeshCreator::default_creating_mesh()
{

    Mesh_criteria criteria(CGAL::parameters::facet_angle   =25.0, 
                           CGAL::parameters::edge_size     =0.025,
                           CGAL::parameters::facet_size    =0.05,
                           CGAL::parameters::facet_distance=0.005,
                           CGAL::parameters::cell_radius_edge_ratio=3,
                           CGAL::parameters::cell_size=0.05);

    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);

    while ( remove_isolated_vertices(c3t3) >0 )
            c3t3.rescan_after_load_of_triangulation();




}
void CGALMeshCreator::create_mesh()
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

void CGALMeshCreator::create_mesh(const double mesh_resolution )
{

    double r = min_sphere.get_bounding_sphere_radius(); 


    const double cell_size = r/mesh_resolution;

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




void CGALMeshCreator::save(std::string OutPath)
{
    std::ofstream  medit_file(OutPath);

    c3t3.output_to_medit(medit_file,false);
    medit_file.close();


}


void CGALMeshCreator::refine_mesh()
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
void CGALMeshCreator::refine_mesh(const double mesh_resolution )
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

void CGALMeshCreator::label_boundary_cells(int btag , int ntag ) // workaround to mark boundary cells of for example lateral ventircles. this allows for easy marking of Facetfunction in FEniCS 
{
  Subdomain_index subdomain_index(btag);
  Subdomain_index subdomain_index_bis(ntag);
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin(subdomain_index);cit != c3t3.cells_in_complex_end(); ++cit)
  {
     for (std::size_t i = 0; i < 4; i++)
     { if (c3t3.subdomain_index(cit->neighbor(i))!=subdomain_index){c3t3.set_subdomain_index(cit->neighbor(i), subdomain_index_bis);}    }      
  }

} 

void CGALMeshCreator::remove_label_cells(int tag) //rename
{
  Subdomain_index subdomain_index(tag);
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin(subdomain_index);cit != c3t3.cells_in_complex_end(); ++cit){ c3t3.remove_from_complex(cit);}
}


//----------------------   Have overloaded functions or just one????????????????
void CGALMeshCreator::add_sharp_border_edges(Polyhedron& polyhedron) // no need to expose ?  
{ 

  Polylines polylinput; 
  typedef boost::property_map<Polyhedron, CGAL::edge_is_feature_t>::type EIF_map;
  EIF_map eif = get(CGAL::edge_is_feature, polyhedron);

  CGAL::Polygon_mesh_processing::detect_sharp_edges(polyhedron,80, eif); // -> threshold ?? 
   for( Polyhedron::Edge_iterator he = polyhedron.edges_begin(); he != polyhedron.edges_end() ; ++he)
   {
      if(he->is_feature_edge() ) 
      {
         Polyline_3 polyline;
         polyline.push_back(he->vertex()->point());
         polyline.push_back(he->opposite()->vertex()->point());     
         polylinput.push_back(polyline);
      }    
  }   
  polylines_to_protect(this->borders, polylinput.begin(),  polylinput.end() ); // borders 
}



void CGALMeshCreator::add_sharp_border_edges(CGALSurface& surface) 
{ 
  Polyhedron polyhedron;
  surface.get_polyhedron(polyhedron);
  add_sharp_border_edges(polyhedron);

}


void CGALMeshCreator::set_features()
{ 
  Polylines polylines;
  polylines_to_protect(polylines, features.begin() , features.end()  );
  set_features(polylines.begin(), polylines.end());
}


void CGALMeshCreator::lloyd(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze )
{CGAL::lloyd_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 

void CGALMeshCreator::odt(double time_limit, int max_iteration_number, double convergence,double freeze_bound, bool do_freeze) 
{CGAL::odt_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } 


#endif
