#include <pybind11/pybind11.h>
/* #include <pybind11/operators.h> */

#include "CGALSurface.h"
#include "CGALMeshCreator.h"
#include "reconstruct_surface.h"


namespace py = pybind11;


PYBIND11_MODULE(brainmesh, m) {

    m.def("reconstruct_surface", &reconstruct_surface);
            /* , */
        /* py::arg("sm_angle") = 20.0, */
        /* py::arg("sm_radius") = 100.0, */
        /* py::arg("sm_distance") = 0.25, */
        /* py::arg("approximation_ratio") = 0.02, */
        /* py::arg("average_spacing_ratio") = 5.0); */

    py::class_<CGALSurface>(m, "BrainSurface")
        .def(py::init<std::string &>())

        /* .def(py::self += py::self) */    // FIXME: Some const trouble

        // Instead of operators
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
        /* .def("points_inside", &CGALSurface::points_inside) */
        /* .def("points_outside", &CGALSurface::points_outside) */

        .def("self_intersections", &CGALSurface::self_intersections)
        .def("num_self_intersections", &CGALSurface::num_self_intersections)

        .def("save", &CGALSurface::save)
        .def("collapse_edges", &CGALSurface::collapse_edges)
        .def("preprocess", &CGALSurface::preprocess)

        /* .def("fair", &CGALSurface::fair) */

        // TODO
        /* .def("insert_surface", &CGALSurface::insert_surface) */  // TODO cpp side
        /* .def("getMesh", &CGALSurface::get_mesh) */       // No need to expose
        /* .def("get_polyhedron", &CGALSurface::get_polyhedron); */ // No need to expose

        // Experimental reconstruction -- creates self-intersections
        .def("reconstruct_surface", &CGALSurface::reconstruct_surface,
                py::arg("sm_angle") = 20,
                py::arg("sm_radius") = 30,
                py::arg("sm_distance") = 0.375)

        .def("num_faces", &CGALSurface::num_faces)
        .def("num_edges", &CGALSurface::num_edges)
        .def("num_vertices", &CGALSurface::num_vertices);

    py::class_<CGALMeshCreator>(m, "BrainMesh")
        .def(py::init<CGALSurface &>())

        // TODO: What to do about theese two? Need more classes?
        /* .def(py::init<std::vector<CGALSurface>, CGAL::Bbox_3, abstract_map>()) */
        /* .def("lipschitz_size_field", &CGALMeshCreator::lipschitz_size_field) */

        .def("set_parameters", &CGALMeshCreator::set_parameters)
        .def("set_parameter", &CGALMeshCreator::set_parameter)
        .def("create_mesh", &CGALMeshCreator::create_mesh)
        .def("default_parameters", &CGALMeshCreator::default_parameters)
        .def("save_mesh", &CGALMeshCreator::save_mesh);
}
