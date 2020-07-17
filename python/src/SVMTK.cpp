#include <pybind11/pybind11.h>
/* #include <pybind11/operators.h> */
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

// // CGAL
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include "Slice.h"
#include "Surface.h"
#include "Domain.h"
#include "surface_mesher.h"

#include "convex_hull.h"



namespace py = pybind11;


class PyAbstractMap : public AbstractMap {
    public:
        using AbstractMap::AbstractMap;  /* Inherit constructors */
};

typedef std::function< double(double, double, double) > Surface_implicit_function;

//using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
//using Pont_3 = Kernel::Point_3;
//using CGAL_point_vector = std::vector< Point_3 >;

// PYBIND11_MAKE_OPAQUE(CGAL_point_vector);

// We want to construct Polyline in Python then let polyline translate into cgal Point in c++

Point_3 Wrapper_point_3(double x, double y, double z)
{
    return Point_3(x, y, z);
}


std::shared_ptr< Surface > convex_hull_wrapper(py::array_t< double > point3_array)
{
    py::buffer_info buffer = point3_array.request();

    // Some sanity checks
    if (buffer.ndim != 2)
        throw std::runtime_error("Expected 2d array");

    if (buffer.size % 3 != 0)
        throw std::runtime_error("Point array size must be shape (N, 3)");

    double *ptr_point3_array = (double *) buffer.ptr;

    auto point_vector = CGAL_point_vector(buffer.size / 3);

    size_t point_vector_counter = 0;
    for (size_t i = 0; i < buffer.size; i += 3)
    {
        double x = ptr_point3_array[i];
        double y = ptr_point3_array[i + 1];
        double z = ptr_point3_array[i + 2];

        point_vector[point_vector_counter++] = Point_3(x, y, z);
    }

    auto surface = convex_hull(point_vector);
    return surface;
}


PYBIND11_MODULE(SVMTK, m) {

    py::class_<Point_3>(m, "Point_3")
       .def(py::init<double, double, double>())
       .def("x", &Point_3::x);

    py::class_<AbstractMap, PyAbstractMap> abstractmap(m, "AbstractMap");

    py::class_<SubdomainMap,AbstractMap>(m, "SubdomainMap")
        .def(py::init<>())
        .def("print",  &SubdomainMap::print)
        .def("add", &SubdomainMap::add);

    py::class_<Slice,std::shared_ptr<Slice>>(m, "Slice")
        .def(py::init<>())
        .def(py::init<Slice&>())
        .def("create_mesh", &Slice::create_mesh)
        .def("simplify", &Slice::simplify)
        .def("keep_component", &Slice::keep_component)
        .def("save", &Slice::save)
        .def("add_constraints",(void (Slice::*)(Slice&,bool)) &Slice::add_constraints);

    py::class_<Surface,std::shared_ptr<Surface>>(m, "Surface")
        .def(py::init<std::string &>())
        .def(py::init<>())

        .def("implicit_surface", &Surface::implicit_surface<Surface_implicit_function>)
        .def("triangulate_hole",&Surface::triangulate_hole)
        .def("clip", &Surface::clip)
        .def("intersection", &Surface::surface_intersection)
        .def("union", &Surface::surface_union)
        .def("difference", &Surface::surface_difference)
        .def("slice", &Surface::mesh_slice)
        .def("span", &Surface::span)

        .def("fill_holes", &Surface::fill_holes)
        .def("triangulate_faces", &Surface::triangulate_faces)
        .def("isotropic_remeshing", &Surface::isotropic_remeshing)
        .def("adjust_boundary", &Surface::adjust_boundary)

        .def("smooth_laplacian", &Surface::smooth_laplacian)
        .def("smooth_taubin", &Surface::smooth_taubin)

        // Either use these two for operator overloading, or return the vertices
        //.def("inside", &Surface::inside)
        //.def("outside", &Surface::outside)
        .def("make_cube", &Surface::make_cube)
        .def("make_cone", &Surface::make_cone)
        .def("make_cylinder", &Surface::make_cylinder)
        .def("make_sphere", &Surface::make_sphere)
        .def("num_self_intersections", &Surface::num_self_intersections)
        .def("collapse_edges", &Surface::collapse_edges)
        .def("save", &Surface::save)
        .def("split_edges", &Surface::split_edges)
        .def("extension", &Surface::cylindric_extension)

        //.def("load", &Surface::load)//
        .def("separate_narrow_gaps", &Surface::seperate_narrow_gaps)
        .def("reconstruct", &Surface::reconstruct)
        .def("convex_hull", &Surface::convex_hull)

        .def("num_faces", &Surface::num_faces)
        .def("num_edges", &Surface::num_edges)
        .def("num_vertices", &Surface::num_vertices);

    py::class_<Domain,std::shared_ptr<Domain>>(m, "Domain")
        .def(py::init<Surface &>())
        .def(py::init<std::vector<Surface>>())
        .def(py::init<std::vector<Surface>, AbstractMap&>())
        .def("set_parameters", &Domain::set_parameters) // std::map<std::string, double>
        .def("set_parameter", &Domain::set_parameter)

        .def("create_mesh", (void (Domain::*)()) &Domain::create_mesh)
        .def("create_mesh", (void (Domain::*)(double)) &Domain::create_mesh)
        .def("refine_mesh", (void (Domain::*)()) &Domain::refine_mesh)
        .def("refine_mesh", (void (Domain::*)(double)) &Domain::refine_mesh)

        .def("get_boundary", &Domain::get_boundary)

        .def("lloyd", &Domain::lloyd)
        .def("odt", &Domain::odt)
        .def("exude", &Domain::exude)
        .def("perturb", &Domain::perturb)

        .def("add_sharp_border_edges", (void (Domain::*)(Surface&,double)) &Domain::add_sharp_border_edges)
        .def("reset_borders", &Domain::reset_borders)
        .def("remove_subdomain", (void (Domain::*)(std::vector<int>)) &Domain::remove_subdomain)
        .def("remove_subdomain", (void (Domain::*)(int)) &Domain::remove_subdomain)
        .def("add_corners",  &Domain::add_corners)
        .def("number_of_cells", &Domain::number_of_cells)
        // add, remove, set
        .def("set_borders", &Domain::set_borders)
        .def("set_features", (void(Domain::*)()) &Domain::set_features)
        .def("set_features", (void (Domain::*)(Surface&)) &Domain::set_features)
        .def("add_feature", &Domain::add_feature)
        .def("save", &Domain::save);

    m.def("separate_surfaces",  (void (*)(Surface&,Surface&,Surface&)) &surface_overlapp<Surface> );
    m.def("separate_surfaces",  (void (*)(Surface&,Surface&)) &surface_overlapp<Surface> );
    m.def("morphological_surface_union", &morphological_surface_union<Surface>);
    m.def("convex_hull", &convex_hull_wrapper);
}
