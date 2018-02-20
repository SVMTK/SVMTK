#ifndef __CGAL_MESH_CREATOR_H


#define __CGAL_MESH_CREATOR_H



//LOCAL
#include "SubdomainMap.h"
/* #include "../CGALSurface/CGALSurface.cpp" */
#include "Polyhedral_vector_to_labeled_function_wrapper.h"



//CGAL

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_3/experimental/Lipschitz_sizing_polyhedron.h>


#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/remove_far_points_in_mesh_3.h>
#include <CGAL/config.h>
#include <CGAL/assertions.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>




//Boost
//C++
#include <map>
#include <list>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

class CGALSurface;

class CGALMeshCreator
{
      public :
           typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
           typedef K::Point_3 Point_3;
           typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
           typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Polyhedral_mesh_domain_3;

           typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, K  > Function_wrapper; //
           typedef Function_wrapper::Function_vector Function_vector; //
           typedef CGAL::Mesh_domain_with_polyline_features_3<CGAL::Labeled_mesh_domain_3<Function_wrapper, K> >  Mesh_domain;

           typedef CGAL::Mesh_3::Lipschitz_sizing<K, Mesh_domain> Lip_sizing;
           typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

           typedef Polyhedral_mesh_domain_3::Corner_index Corner_index;
           typedef Polyhedral_mesh_domain_3::Curve_segment_index Curve_segment_index;

           typedef CGAL::Mesh_complex_3_in_triangulation_3< Tr, Corner_index, Curve_segment_index> C3t3;

           typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;


           typedef std::map< std::string ,double> Parameters;

            CGALMeshCreator(std::vector<CGALSurface> surfaces, CGAL::Bbox_3 bbox_3 , AbstractMap &map);

            CGALMeshCreator(std::vector<CGALSurface> surfaces, AbstractMap& map);

            CGALMeshCreator(CGALSurface& surface,  CGAL::Bbox_3 bbox_3 );

            CGALMeshCreator(CGALSurface& surface);

            ~CGALMeshCreator() {}

            void add_polylines(const std::string filename);

            template< typename InputIterator>
            void add_polylines(InputIterator begin, InputIterator end); //wrapper double to points

            void lipschitz_size_field(int subdomain_id, int k,double min_size,double max_size);

            void set_parameters(Parameters new_parameters);

            void set_parameter(std::string key , double value );

            void create_mesh();

            void create_mesh( int initial_points);

            void default_parameters()
            {
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
//          void remove_isolated_vertices(); stand alone function ??
            void save_mesh(const std::string OutPath);
      private :
            //const Labeled_mesh_domain_3& r_domain_;
            std::unique_ptr<Mesh_domain> domain_ptr;
            std::unique_ptr<Lip_sizing>  lip_sizing_ptr;
            Parameters parameters;
            C3t3 c3t3;
              
};




#endif
