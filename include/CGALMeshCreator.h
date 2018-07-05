#ifndef __CGAL_MESH_CREATOR_H
#define __CGAL_MESH_CREATOR_H

#define BOOST_PARAMETER_MAX_ARITY 12

#include "CGALSurface.h"

/* //LOCAL */
/* #include "SubdomainMap.h" */
/* #include "../CGALSurface/CGALSurface.cpp" */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>     // Kernel
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_polyhedron_3.h>

#include "Polyhedral_vector_to_labeled_function_wrapper.h"

#include <CGAL/Labeled_mesh_domain_3.h>     //
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_3/experimental/Lipschitz_sizing_polyhedron.h>

// Mesh criteria
#include <CGAL/Mesh_criteria_3.h>

// Mesh_complex_3_in_triangulation_3
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/make_mesh_3.h>

/* #include <CGAL/Mesh_criteria_3.h> */
/* #include <CGAL/Mesh_cell_base_3.h> */
/* #include <CGAL/Mesh_vertex_base_3.h> */
/* #include <CGAL/remove_far_points_in_mesh_3.h> */
/* #include <CGAL/config.h> */
/* #include <CGAL/assertions.h> */

/* #include <CGAL/IO/File_medit.h> */
/* #include <CGAL/make_mesh_3.h> */
/* #include <CGAL/refine_mesh_3.h> */
/* #include <CGAL/IO/Polyhedron_iostream.h> */




/* typedef CGAL::Mesh_3::Lipschitz_sizing<K, Mesh_domain> Lip_sizing; */





class CGALMeshCreator {
    typedef std::map<std::string, double> Parameters;

    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef Kernel::Point_3 Point_3;
    typedef CGAL::Mesh_polyhedron_3<Kernel>::type Polyhedron;
    typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel> Polyhedral_mesh_domain_3;
    typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, Kernel> Function_wrapper; //
    typedef Function_wrapper::Function_vector Function_vector; //
    typedef CGAL::Mesh_domain_with_polyline_features_3<CGAL::Labeled_mesh_domain_3<Function_wrapper, Kernel> >  Mesh_domain;

    typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
    typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

    typedef Polyhedral_mesh_domain_3::Corner_index Corner_index;
    typedef Polyhedral_mesh_domain_3::Curve_segment_index Curve_segment_index;
    typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, Corner_index, Curve_segment_index> C3t3;

    Parameters parameters;
        std::unique_ptr<Mesh_domain> domain_ptr;
        /* std::unique_ptr<Lip_sizing>  lip_sizing_ptr; */
        C3t3 c3t3;

    public:
        CGALMeshCreator(CGALSurface &);

        void set_parameters(const Parameters &);

        void set_parameter(const std::string, const double);

        void create_mesh(const int initial_points = 0);

        void default_parameters();
};


CGALMeshCreator::CGALMeshCreator(CGALSurface &surface) {
    default_parameters();

    Function_vector funcvec;
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);

    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
    funcvec.push_back(polyhedral_domain);
    Function_wrapper wrapper(funcvec);
    Mesh_domain domain(wrapper, wrapper.bbox());
    domain_ptr = std::unique_ptr<Mesh_domain>(new Mesh_domain(wrapper, wrapper.bbox()));
}


void CGALMeshCreator::default_parameters() {
    parameters["mesh_resolution"] = 64.0;
    parameters["perturb_optimize"] = 0.0;
    parameters["exude_optimize"] = 0.0;
    parameters["lloyd_optimize"] = 0.0;
    parameters["odt_optimize"] = 0.0;
    parameters["edge_size"] = 0.25;
    parameters["facet_angle"] = 25.0;
    parameters["facet_size"] = 0.1;
    parameters["facet_distance"] = 0.1;
    parameters["cell_radius_edge_ratio"] = 3.0;
    parameters["cell_size"] = 0.1;
    parameters["detect_sharp_features"] = 1.0;
    parameters["feature_threshold"]= 70.;
}


void CGALMeshCreator::set_parameters(const Parameters &new_parameters) {
    for (const auto &kv: new_parameters) {
        parameters[kv.first] = static_cast<double>(kv.second);
        std::cout << kv.first << ": " << kv.second << std::endl;
    }
}


void CGALMeshCreator::set_parameter(const std::string key, const double value) {
    parameters[key] = value;
}


void CGALMeshCreator::create_mesh(const int initial_points) {
    Mesh_criteria criteria(
            CGAL::parameters::edge_size=parameters["edge_size"],
            CGAL::parameters::facet_angle=parameters["facet_angle"],
            CGAL::parameters::facet_size=parameters["facet_size"],
            CGAL::parameters::facet_distance=parameters["facet_distance"],
            CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
            CGAL::parameters::cell_size=parameters["cell_size"]);

    /* CGAL::internal::Mesh_3::init_c3t3(c3t3, *domain_ptr.get(), criteria, initial_points); */
    C3t3 c3t3= CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);
    /* refine_mesh_3(c3t3, *domain_ptr.get(), criteria, CGAL::parameters::no_reset_c3t3()); */

    /* if(parameters["odt_optimize"] ) { */
    /* } */

    /* if(parameters["perturb_optimize"]) { */
    /* } */
    /* if(parameters["lloyd_optimize"]) { */
    /* } */
    /* if( parameters["exude_optimize"] ) { */
    /* } */
}

#endif
