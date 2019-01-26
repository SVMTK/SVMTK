#include <pybind11/pybind11.h>
/* #include <pybind11/operators.h> */
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include "CGALSurface.h"
#include "CGALMeshCreator.h"
#include "Neuron.h"
//#include "reconstruct_surface.h" 


namespace py = pybind11;


class PyAbstractMap : public AbstractMap{
public:
       using AbstractMap::AbstractMap; /* Inherit constructors */

};

PYBIND11_MODULE(brainmesh, m) {

    /* m.def("reconstruct_surface", &reconstruct_surface); */
            /* , */
        /* py::arg("sm_angle") = 20.0, */
        /* py::arg("sm_radius") = 100.0, */
        /* py::arg("sm_distance") = 0.25, */
        /* py::arg("approximation_ratio") = 0.02, */
        /* py::arg("average_spacing_ratio") = 5.0); */



      // m.def("reconstruct_surface", &reconstruct_surface );


    py::class_<AbstractMap,PyAbstractMap> abstractmap(m,"AbstractMap");

    py::class_<SubdomainMap>(m, "SubdomainMap")
        .def(py::init<>())
        .def("__repr__",  &SubdomainMap::print)
        .def("add", &SubdomainMap::add);

    py::class_<CGALSurface>(m, "BrainSurface")
        .def(py::init<std::string &>())
        .def(py::init<>())
        /* .def(py::self += py::self) */    // FIXME: Some const trouble



        .def("intersection", &CGALSurface::surface_intersection)
        .def("union", &CGALSurface::surface_union)
        .def("difference", &CGALSurface::surface_difference)

        .def("fill_holes", &CGALSurface::fill_holes)
        .def("triangulate_faces", &CGALSurface::triangulate_faces)
        .def("stitch_borders", &CGALSurface::stitch_borders)
        .def("isotropic_remeshing", &CGALSurface::isotropic_remeshing)
        .def("adjust_boundary", &CGALSurface::adjust_boundary)

        .def("smooth_laplacian", &CGALSurface::smooth_laplacian)
        .def("smooth_taubin", &CGALSurface::smooth_taubin)

        // Either use these two for operator overloading, or return the vertices
        .def("points_inside", &CGALSurface::points_inside)
        .def("points_outside", &CGALSurface::points_outside)
        .def("make_cube", &CGALSurface::make_cube)
	.def("make_cone", &CGALSurface::make_cone)
	.def("make_cylinder", &CGALSurface::make_cylinder)
        .def("make_sphere", &CGALSurface::make_sphere)

        .def("self_intersections", &CGALSurface::self_intersections)
        .def("num_self_intersections", &CGALSurface::num_self_intersections) 
        .def("collapse_edges", &CGALSurface::collapse_edges)
        .def("save", &CGALSurface::save)

        .def("fair", &CGALSurface::fair)
        //.def("load", &CGALSurface::load)//
        .def("fix_close_junctures", &CGALSurface::fix_close_junctures)
        .def("reconstruct_surface", &CGALSurface::reconstruct_surface)

        /* // Experimental reconstruction -- creates self-intersections */
        /* .def("reconstruct_surface", &CGALSurface::reconstruct_surface, */
        /*         py::arg("sm_angle") = 20, */
        /*         py::arg("sm_radius") = 30, */
        /*         py::arg("sm_distance") = 0.375) */

        .def("num_faces", &CGALSurface::num_faces)
        .def("num_edges", &CGALSurface::num_edges)
        .def("num_vertices", &CGALSurface::num_vertices);

    py::class_<CGALMeshCreator>(m, "BrainMesh")
        .def(py::init<CGALSurface &>())
        .def(py::init<std::vector<CGALSurface>>())
        .def(py::init<std::vector<CGALSurface>, AbstractMap&>())

        .def("create_mesh", (void (CGALMeshCreator::*)()) &CGALMeshCreator::create_mesh) 
        .def("create_mesh", (void (CGALMeshCreator::*)(int)) &CGALMeshCreator::create_mesh)

        .def("lloyd", &CGALMeshCreator::lloyd)
        .def("odt", &CGALMeshCreator::odt)
        .def("excude", &CGALMeshCreator::excude)
        .def("perturb", &CGALMeshCreator::perturb)

        .def("add_sharp_border_edges", (void (CGALMeshCreator::*)(CGALSurface&)) &CGALMeshCreator::add_sharp_border_edges)
        .def("refine_mesh", &CGALMeshCreator::refine_mesh)


    /*     // TODO: What to do about theese two? Need more classes? */
    /*     /1*  *1/ */
    /*     /1* .def("lipschitz_size_field", &CGALMeshCreator::lipschitz_size_field) Make subclass with lipschits *1/ */ 

        .def("set_parameters", &CGALMeshCreator::set_parameters) // std::map<std::string, double>
        .def("set_parameter", &CGALMeshCreator::set_parameter)


        .def("set_borders", &CGALMeshCreator::set_borders)

        .def("set_features", (void(CGALMeshCreator::*)(CGALMeshCreator::Polylines&)) &CGALMeshCreator::set_features) 
        .def("save_mesh", &CGALMeshCreator::save_mesh); 


    py::class_<Neuron,CGALSurface>(m,"Neuron")
        .def(py::init<std::string>())
        .def("surface_mesh", &Neuron::surface_mesh );




       m.def("surface_overlapp", &surface_overlapp<CGALSurface> );


}
