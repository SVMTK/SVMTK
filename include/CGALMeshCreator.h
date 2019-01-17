#ifndef __CGAL_MESH_CREATOR_H


#define __CGAL_MESH_CREATOR_H


#include "CGALSurface.h" 
#include "SubdomainMap.h" 
#include "Polyhedral_vector_to_labeled_function_wrapper.h"

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // CGALSURFACE  CGALSURFACE ? 
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_3/experimental/Lipschitz_sizing_polyhedron.h>
//#include <CGAL/Mesh_3/Detect_features_in_polyhedra.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/remove_far_points_in_mesh_3.h>
#include <CGAL/config.h>
#include <CGAL/assertions.h>

//#include <CGAL/IO/File_medit.h>
//#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
//#include <CGAL/IO/Polyhedron_iostream.h>

// 	Mesh_facet_topology parameters::facet_topology  	= CGAL::FACET_VERTICES_ON_SURFACE, 
//Boost
//C++


//#include <map>

#include <list>
//#include <string>
//#include <vector>
//#include <iostream>
#include <fstream>
#include <memory>



// CURRENT PLAN :
// TODO: Consider void remove_isolated_vertices() as stand alone function or class function
// TODO: Clean up header files
// TODO: Consider adding MeshCriteria as class member 
// TODO: Fix Lip_sizing field 



class CGALSurface;

class CGALMeshCreator {
    public :
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef K::Point_3 Point_3;
        typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron; // ???

        // temporary -- Until I can get Lars code to work
        /* typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain; */
        // End temporary


        typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron> Polyhedral_mesh_domain_3; // ?? vs CGAL::Polyhedral_mesh_domain_with_features_3<K> Polyhedral_mesh_domain_3;

        // FIXME: There is a problem with the function wrapper
        typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, K  > Function_wrapper; //
        typedef Function_wrapper::Function_vector Function_vector; //
        typedef CGAL::Labeled_mesh_domain_3< K> Labeled_Mesh_Domain;
        typedef CGAL::Mesh_domain_with_polyline_features_3<Labeled_Mesh_Domain> Mesh_domain; // labeled mesh_domain function wrapper 

        typedef CGAL::Mesh_3::Lipschitz_sizing<K, Mesh_domain> Lip_sizing;
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


        CGALMeshCreator(CGALSurface& surface);
        CGALMeshCreator(std::vector<CGALSurface> surfaces, CGAL::Bbox_3 bbox_3 , AbstractMap &map);
        CGALMeshCreator(std::vector<CGALSurface> surfaces, AbstractMap& map);
        CGALMeshCreator(std::vector<CGALSurface> surfaces);
        CGALMeshCreator(CGALSurface& surface,  CGAL::Bbox_3 bbox_3 );


        ~CGALMeshCreator() {}

            void lipschitz_size_field(int subdomain_id, int k,double min_size,double max_size);

        void set_parameters(Parameters new_parameters);

        void set_parameter(std::string key, double value);

        void create_mesh();

            void create_mesh( int initial_points);

        void default_parameters() {
            parameters["mesh_resolution"]=64.0;
            parameters["perturb_optimize"] =0.0;
            parameters["exude_optimize"] = 0.0;
            parameters["lloyd_optimize"] =  0.0;
            parameters["odt_optimize"]   =  0.0;
            parameters["edge_size"]      =   0.25;
            parameters["facet_angle"]    = 25.0;
            parameters["facet_size"]     = 0.1;
            parameters["facet_distance"] =  0.1;
            parameters["cell_radius_edge_ratio"] = 3.0;
            parameters["cell_size"] = 0.1;
            parameters["detect_sharp_features"] = 1.0;
            parameters["feature_threshold"]= 70.;
        }

        void save_mesh(std::string OutPath);

            void add_sharp_edges(Polyhedron& polyhedron) ;

            void refine_mesh();

            Polylines& get_polylines() {return pedges; }

            void add_polyline(Polyline_3 polyline) { pedges.push_back(polyline);}   

            void reset_polyline(){pedges.clear();}

            void set_polylines(){domain_ptr.get()->add_features(get_polylines().begin(), get_polylines().end());} 

            template< typename InputIterator> // Polylines::iterator
            void set_polylines(InputIterator begin, InputIterator end){domain_ptr->add_features(begin, end);}

            template< typename InputIterator>
            void insert_edges(InputIterator begin, InputIterator end); 

            template < typename InputIterator>
            void insert_edge(InputIterator begin, InputIterator end);

            /* // Built due to private mesh, better option ? */ 
            /* void lloyd(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true) */
            /* {CGAL::lloyd_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } */ 

            /* void excude( double time_limit = 0, double sliver_bound = 0 ){ CGAL::exude_mesh_3(c3t3, sliver_bound=sliver_bound, time_limit=time_limit);} */

            /* void perturb( double time_limit=0, double sliver_bound=0){CGAL::perturb_mesh_3 ( c3t3, *domain_ptr.get(), time_limit=time_limit, sliver_bound=sliver_bound) ;} */

            /* void odt(double time_limit= 0, int max_iteration_number = 0, double convergence = 0.02,double freeze_bound = 0.01, bool do_freeze = true) */ 
            /* {CGAL::odt_optimize_mesh_3(c3t3, *domain_ptr.get(), time_limit=time_limit, max_iteration_number=max_iteration_number,convergence=convergence, freeze_bound  = freeze_bound, do_freeze = do_freeze); } */ 
            //

            void label_boundary_cells(int btag, int ntag); // Work around 

            void remove_label_cells(int tag);    // tags 

            //void dolfin_mesh(); require linking to dolfin 

    private :
        std::unique_ptr<Mesh_domain> domain_ptr;
        std::unique_ptr<Lip_sizing>  lip_sizing_ptr;
        Polylines  pedges; 
        Parameters parameters;
        C3t3 c3t3;
};




#endif
