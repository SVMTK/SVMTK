#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

//#include "CGALMeshCreator.h"
#include "CGALSurface.h"
#include "implicit_functions.h"
#include "Polyhedral_vector_to_labeled_function_wrapper.h"
#include "read_polygons_STL.h"
#include "read_polylines.h"
#include "remove_isolated_vertices.h"
#include "SubdomainMap.h"

#include "CGALMeshCreator.h"


namespace py = pybind11;


//Trampoline class to tell pybind how to deal with virutal functions
/* class pyAbstractMap : public AbstractMap { */
/*     public: */ 
/*         using AbstractMap::AbstractMap;     // Inherit constructors */
/*         int index(boost::dynamic_bitset<> bits) override { */
/*             PYBIND11_OVERLOAD_PURE( */
/*                int, */
/*                AbstractMap, */
/*                index, */
/*                bits */
/*             ); */
/*         } */
/* }; */



PYBIND11_MODULE(brainmesh, m) {
    py::class_<CGALSurface>(m, "BrainSurface")
        .def(py::init<std::string &>())
        /* .def(py::init<const Implicit_function, double, double, double, double>()) */
        /* .def(py::self += py::self) */ 
        .def("fill_holes", &CGALSurface::fill_holes)
        .def("triangulate_faces", &CGALSurface::triangulate_faces)
        .def("stitch_borders", &CGALSurface::stitch_borders)
        .def("insert_surface", &CGALSurface::insert_surface)
        .def("isotropic_remeshing", &CGALSurface::isotropic_remeshing)
        .def("adjust_boundary", &CGALSurface::adjust_boundary)
        .def("smooth_laplacian", &CGALSurface::smooth_laplacian)
        /* .def("adjusting_boundary_region", &CGALSurface::adjusting_boundary_region) */
        /* .def("smooth_laplacian_region", &CGALSurface::smooth_laplacian_region) */
        .def("points_inside", &CGALSurface::points_inside)
        .def("points_outside", &CGALSurface::points_outside)
        .def("getMesh", &CGALSurface::get_mesh)
        .def("self_intersections", &CGALSurface::self_intersections)
        /* .def("get_polyhedron", &CGALSurface::get_polyhedron) */
        .def("save", &CGALSurface::save)
        .def("preprocess", &CGALSurface::preprocess)
        .def("fair", &CGALSurface::fair);

    /* py::class_<AbstractMap, pyAbstractMap> abstract_map(m, "abstract_map"); */
    /* abstract_map */
    /*     .def(py::init<>()) */
    /*     .def("index", &AbstractMap::index); */

    py::class_<CGALMeshCreator>(m, "BrainMesh")
        /* .def(py::init<std::vector<CGALSurface>, CGAL::Bbox_3, abstract_map>()) */
        /* .def(py::init<CGALSurface &>()) */
        /* .def(py::init<CGALSurface &>()) */
        .def(py::init<CGALSurface &>())
        .def("lipschitz_size_field", &CGALMeshCreator::lipschitz_size_field)
        .def("set_parameters", &CGALMeshCreator::set_parameters)
        .def("set_parameter", &CGALMeshCreator::set_parameter)
        .def("create_mesh", &CGALMeshCreator::create_mesh)
        /* .def("create_mesh", &CGALMeshCreator::create_mesh) */    // Overloaded function. Look up
        .def("default_parameters", &CGALMeshCreator::default_parameters)
        .def("save_mesh", &CGALMeshCreator::save_mesh);
}
