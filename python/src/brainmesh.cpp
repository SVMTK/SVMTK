#include <pybind11/pybind11.h>
/* #include <pybind11/operators.h> */

#include "CGALSurface.h"
#include "CGALMeshCreator.h"
#include "reconstruct_surface.h"


namespace py = pybind11;


class PyAbstractMap : public AbstractMap{
    public:
        using AbstractMap::AbstractMap; /* Inherit constructors */
};


typedef std::function<double(double,double,double)> Surface_implicit_function;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;


// We want to construct Polyline in Python then let polyline translate into cgal Point in c++

Point_3 Wrapper_point_3(const double x, const double y, const double z) {
    return Point_3(x, y, z);
}


PYBIND11_MODULE(brainmesh, m) {

    m.def("reconstruct_surface", &reconstruct_surface,
        py::arg("input_surface"),
        py::arg("sm_angle") = 20.0,
        py::arg("sm_radius") = 100.0,
        py::arg("sm_distance") = 0.25,
        py::arg("approximation_ratio") = 0.02,
        py::arg("average_spacing_ratio") = 5.0);

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
        .def("points_inside", &CGALSurface::points_inside)
        .def("points_outside", &CGALSurface::points_outside)

        .def("self_intersections", &CGALSurface::self_intersections)
        .def("num_self_intersections", &CGALSurface::num_self_intersections)

        .def("save", &CGALSurface::save)
        .def("collapse_edges", &CGALSurface::collapse_edges)
        .def("preprocess", &CGALSurface::preprocess)

        .def("fair", &CGALSurface::fair)

        .def("make_cone", &CGALSurface::make_cone)
        .def("make_cylinder", &CGALSurface::make_cylinder)
        .def("make_cube", &CGALSurface::make_cube)

        .def("fix_close_junctures", &CGALSurface::fix_close_junctures)

        // TODO
        /* .def("insert_surface", &CGALSurface::insert_surface) */  // TODO cpp side
        /* .def("getMesh", &CGALSurface::get_mesh) */       // No need to expose
        /* .def("get_polyhedron", &CGALSurface::get_polyhedron); */ // No need to expose

        .def("num_faces", &CGALSurface::num_faces)
        .def("num_edges", &CGALSurface::num_edges)
        .def("num_vertices", &CGALSurface::num_vertices);

    py::class_<CGALMeshCreator>(m, "BrainMesh")
        .def(py::init<CGALSurface &>())

        .def("create_mesh", &CGALMeshCreator::create_mesh)
        .def("refine_mesh", &CGALMeshCreator::refine_mesh)
        .def("create_mesh", (void (CGALMeshCreator::*)()) &CGALMeshCreator::create_mesh)
        // .def("create_mesh", (void (CGALMeshCreator::*)(int)) &CGALMeshCreator::create_mesh)
        .def("create_mesh", (void (CGALMeshCreator::*)(double)) &CGALMeshCreator::create_mesh)
        .def("default_creating_mesh", &CGALMeshCreator::default_creating_mesh)

        .def("lloyd", &CGALMeshCreator::lloyd)
        .def("odt", &CGALMeshCreator::odt)
        .def("excude", &CGALMeshCreator::excude)
        .def("perturb", &CGALMeshCreator::perturb)

        .def("add_sharp_border_edges", (void (CGALMeshCreator::*)(CGALSurface&)) &CGALMeshCreator::add_sharp_border_edges)

        /* .def("refine_mesh", &CGALMeshCreator::refine_mesh) */
        .def("reset_borders", &CGALMeshCreator::reset_borders)

    /*     // TODO: What to do about theese two? Need more classes? */
    /*     /1* .def(py::init<std::vector<CGALSurface>, CGAL::Bbox_3, abstract_map>()) *1/ */
    /*     /1* .def("lipschitz_size_field", &CGALMeshCreator::lipschitz_size_field) *1/ */

    /*     .def("set_parameters", &CGALMeshCreator::set_parameters) */
    /*     .def("set_parameter", &CGALMeshCreator::set_parameter) */
    /*     .def("create_mesh", &CGALMeshCreator::create_mesh) */
    /*     .def("default_parameters", &CGALMeshCreator::default_parameters) */
        .def("save_mesh", &CGALMeshCreator::save_mesh);
}
