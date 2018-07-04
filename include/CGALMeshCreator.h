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




/* typedef CGAL::Mesh_3::Lipschitz_sizing<K, Mesh_domain> Lip_sizing; */
/* typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr; */

/* typedef Polyhedral_mesh_domain_3::Corner_index Corner_index; */
/* typedef Polyhedral_mesh_domain_3::Curve_segment_index Curve_segment_index; */

/* typedef CGAL::Mesh_complex_3_in_triangulation_3< Tr, Corner_index, Curve_segment_index> C3t3; */

/* typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria; */


class CGALMeshCreator {
    typedef std::map<std::string, double> Parameters;

    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    typedef Kernel::Point_3 Point_3;
    typedef CGAL::Mesh_polyhedron_3<Kernel>::type Polyhedron;
    typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel> Polyhedral_mesh_domain_3;
    typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, Kernel> Function_wrapper; //
    typedef Function_wrapper::Function_vector Function_vector; //
    typedef CGAL::Mesh_domain_with_polyline_features_3<CGAL::Labeled_mesh_domain_3<Function_wrapper, Kernel> >  Mesh_domain;

    Parameters parameters;
        std::unique_ptr<Mesh_domain> domain_ptr;
        /* std::unique_ptr<Lip_sizing>  lip_sizing_ptr; */
        /* C3t3 c3t3; */

    public:
        CGALMeshCreator(CGALSurface &);

        /* void set_parameters(Parameters new_parameters) {}; */

        /* void set_parameter(std::string key, double value) {}; */

        /* void create_mesh() {}; */

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

#endif
