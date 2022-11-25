#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>
#include <pybind11/chrono.h>

#include "docstrings.h"
#include "Surface.h"
#include "Domain.h"
#include "Slice.h"
#include <CGAL/IO/read_ply_points.h>

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

class PyAbstractMap : public AbstractMap
{
public:
    using AbstractMap::AbstractMap; /* Inherit constructors */
};

typedef std::function<double(double, double, double)> Surface_implicit_function;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Vector_3 Vector_3;

/* -- CGAL wrapping */
Vector_3 Wrapper_vector_3(double x, double y, double z)
{
    return Vector_3(x, y, z);
}
Point_3 Wrapper_point_3(double x, double y, double z)
{
    return Point_3(x, y, z);
}
Point_2 Wrapper_point_2(double x, double y)
{
    return Point_2(x, y);
}
Plane_3 Wrapper_plane_3(double x1, double x2, double x3, double x4)
{
    return Plane_3(x1, x2, x3, x4);
}
Plane_3 Wrapper_plane_3(Point_3 p1, Vector_3 v1)
{
    return Plane_3(p1, v1);
}
Plane_3 Wrapper_plane_3(Point_3 p1, Point_3 p2, Point_3 p3)
{
    return Plane_3(p1, p2, p3);
}
std::vector<Point_3> Wrapper_load_points(std::string filename)
{
    std::vector<Point_3> points;
    std::ifstream in(filename);
    if (!CGAL::IO::read_PLY(in, std::back_inserter(points)))
        throw InvalidArgumentError("Can't open file.");
    return points;
}

std::shared_ptr<Surface> Wrapper_convex_hull(py::array_t<double> point3_array)
{
    py::buffer_info buffer = point3_array.request();

    if (buffer.ndim != 2)
        throw std::runtime_error("Expected 2d array");

    if (buffer.size % 3 != 0)
        throw std::runtime_error("Point array size must be shape (N, 3)");

    double *ptr_point3_array = (double *)buffer.ptr;

    auto point_vector = std::vector<Surface::Point_3>(buffer.size / 3);

    size_t point_vector_counter = 0;
    for (int i = 0; i < buffer.size; i += 3)
    {
        double x = ptr_point3_array[i];
        double y = ptr_point3_array[i + 1];
        double z = ptr_point3_array[i + 2];

        point_vector[point_vector_counter++] = Point_3(x, y, z);
    }

    auto surface = convex_hull<Surface, Point_3>(point_vector);
    return surface;
}

PYBIND11_MODULE(SVMTK, m)
{
    m.doc() = "Surface Volume Meshing Toolkit";

    static py::exception<PreconditionError> ex1(m, "PreconditionError");
    static py::exception<EmptyMeshError> ex2(m, "EmptyMeshError");
    static py::exception<InvalidArgumentError> ex3(m, "InvalidArgumentError");
    static py::exception<AlgorithmError> ex4(m, "AlgorithmError");

    py::register_exception_translator([](std::exception_ptr p)
                                      {
                                          try
                                          {
                                              if (p)
                                                  std::rethrow_exception(p);
                                          }
                                          catch (const PreconditionError &e)
                                          {
                                              ex1(e.what());
                                          }
                                          catch (const InvalidArgumentError &e)
                                          {
                                              ex3(e.what());
                                          }
                                          catch (const EmptyMeshError &e)
                                          {
                                              ex2(e.what());
                                          }
                                          catch (const AlgorithmError &e)
                                          {
                                              ex4(e.what());
                                          } });

    py::class_<Vector_3, std::shared_ptr<Vector_3>>(m, "Vector_3", DOC(Vector3))
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"), DOC(Vector3, Vector3))
        .def(py::init<Point_3, Point_3>(), py::arg("source"), py::arg("target"), DOC(Vector3, Vector3, 2))
        .def("__repr__", [](Vector_3 const &self)
             {
                 std::ostringstream os;
                 os << "(" << self.x() << ", " << self.y() << ", " << self.z() << ")";
                 return os.str(); })
        .def(float() * py::self)
        .def(py::self * float())
        .def(py::self + py::self)
        .def("squared_length", &Vector_3::squared_length, "Returns the squared length of the Vector")
        .def("x", &Vector_3::x, "Returns x direction.")
        .def("y", &Vector_3::y, "Returns y direction.")
        .def("z", &Vector_3::z, "Returns z direction.");

    py::class_<Plane_3, std::shared_ptr<Plane_3>>(m, "Plane_3", DOC(Plane3))
        .def(py::init<double, double, double, double>(), py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"), DOC(Plane3, Plane3))
        .def(py::init<Point_3, Vector_3>(), py::arg("point"), py::arg("vector"), DOC(Plane3, Plane3, 2))
        .def("__repr__", [](Plane_3 const &self)
             {
                 std::ostringstream os;
                 os << "(" << self.a() << ", " << self.b() << ", " << self.c() << "," << self.d() << ")";
                 return os.str(); })
        .def("a", &Plane_3::a, "Returns first coefficient of plane equation.")
        .def("b", &Plane_3::b, "Returns second coefficient of plane equation.")
        .def("c", &Plane_3::c, "Returns third coefficient of plane equation.")
        .def("d", &Plane_3::d, "Returns fourth coefficient of plane equation.");

    py::class_<Point_3, std::shared_ptr<Point_3>>(m, "Point_3", DOC(Point3))
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"), DOC(Point3, Point3))
        .def("__repr__", [](Point_3 const &self)
             {
           std::ostringstream os;
           os << "(" << self.x() <<", "<< self.y()<<", "<< self.z() <<")";
           return os.str(); })
        .def(
            "__eq__", [](Point_3 const &self, Point_3 const &other)
            { return (self.x() == other.x() and self.y() == other.y() and self.z() == other.z()); },
            py::is_operator())
        .def(
            "__ne__", [](Point_3 const &self, Point_3 const &other)
            { return (self.x() != other.x() and self.y() != other.y() and self.z() != other.z()); },
            py::is_operator())
        .def(py::self + Vector_3())
        .def(py::self += Vector_3())
        .def(py::self - Vector_3())
        .def("x", &Point_3::x, "Returns x coordinate.")
        .def("y", &Point_3::y, "Returns y coordinate.")
        .def("z", &Point_3::z, "Returns z coordinate.");

    py::class_<Point_2, std::shared_ptr<Point_2>>(m, "Point_2", DOC(Point2))
        .def(py::init<double, double>(), DOC(Point2, Point2))
        .def("__repr__", [](Point_2 const &self)
             {
                 std::ostringstream os;
                 os << "(" << self.x() << ", " << self.y() << ")";
                 return os.str(); })
        .def("x", &Point_2::x, "Returns x coordinate.")
        .def("y", &Point_2::y, "Returns y coordinate.");

    py::class_<AbstractMap, PyAbstractMap, std::shared_ptr<AbstractMap>>(m, "AbstractMap");

    py::class_<SubdomainMap, AbstractMap, std::shared_ptr<SubdomainMap>>(m, "SubdomainMap", DOC(SubdomainMap))
        .def(py::init<int>(), py::arg("num_surfaces") = 0, DOC(SubdomainMap, SubdomainMap))
        .def("__repr__", [](SubdomainMap &self) // TODO TEST
             { return self.print(); })
        .def("print", &SubdomainMap::print, DOC(SubdomainMap, print))
        .def("set_number_of_surfaces", &SubdomainMap::set_number_of_surfaces, DOC(SubdomainMap, set_number_of_surfaces))
        .def("add_interface", &SubdomainMap::add_interface, DOC(SubdomainMap, add_interface))
        .def("get_interfaces", &SubdomainMap::get_interfaces, DOC(SubdomainMap, get_interfaces))
        .def("erase", &SubdomainMap::erase, DOC(SubdomainMap, erase))
        .def("get_map", &SubdomainMap::get_map, DOC(SubdomainMap, get_map))
        .def("get_tags", &SubdomainMap::get_tags, DOC(SubdomainMap, get_tags))
        .def("add", &SubdomainMap::add, DOC(SubdomainMap, add));

    py::class_<Slice, std::shared_ptr<Slice>>(m, "Slice", DOC(Slice))
        .def(py::init<>(), DOC(Slice, Slice))
        .def(py::init<Plane_3>(), DOC(Slice, Slice, 2))
        .def(py::init<Point_3, Vector_3>(), DOC(Slice, Slice, 3))
        .def(py::init<double, double, double, double>(), DOC(Slice, Slice, 4))

        .def("create_mesh", py::overload_cast<double>(&Slice::create_mesh), DOC(Slice, create_mesh))
        .def("create_mesh", py::overload_cast<double, double>(&Slice::create_mesh), DOC(Slice, create_mesh, 2))
        .def("simplify", &Slice::simplify, DOC(Slice, simplify))
        .def("save", &Slice::save, DOC(Slice, save))
        .def("slice_surfaces", &Slice::slice_surfaces<Surface>, DOC(Slice, slice_surfaces))
        .def("export_as_surface", &Slice::export_as_surface<Surface>, DOC(Slice, export_as_surface))
        .def("add_surface_domains", py::overload_cast<std::vector<Surface>, AbstractMap &>(&Slice::add_surface_domains<Surface>), DOC(Slice, add_surface_domains))
        .def("add_surface_domains", py::overload_cast<std::vector<Surface>>(&Slice::add_surface_domains<Surface>), DOC(Slice, add_surface_domains, 2))
        .def("number_of_constraints", &Slice::number_of_constraints, DOC(Slice, number_of_constraints))
        .def("number_of_subdomains", &Slice::number_of_subdomains, DOC(Slice, number_of_subdomains))
        .def("number_of_faces", &Slice::number_of_faces, DOC(Slice, number_of_faces))
        .def("connected_components", &Slice::connected_components, DOC(Slice, connected_components))
        .def("keep_largest_connected_component", &Slice::keep_largest_connected_component, DOC(Slice, keep_largest_connected_component))
        .def("remove_subdomain", py::overload_cast<int>(&Slice::remove_subdomain), DOC(Slice, remove_subdomain))
        .def("remove_subdomain", py::overload_cast<std::vector<int>>(&Slice::remove_subdomain), DOC(Slice, remove_subdomain, 2))
        .def("get_constraints", &Slice::get_constraints, DOC(Slice, get_constraints))
        .def("add_constraint", &Slice::add_constraint, DOC(Slice, add_constraint))
        .def("add_constraints", py::overload_cast<Slice &>(&Slice::add_constraints), DOC(Slice, get_constraints));
    // clear constraints

    py::class_<Surface, std::shared_ptr<Surface>>(m, "Surface", DOC(Surface))
        .def(py::init<>(), DOC(Surface, Surface))
        .def(py::init<std::string &>(), py::arg("filename"), DOC(Surface, Surface, 2))
        .def(py::init<Surface &>(), py::arg("surf"), DOC(Surface, Surface, 3))
        .def("__copy__", [](const Surface &self)
             { return Surface(self); })
        .def("copy", [](const Surface &self)
             { return Surface(self); })
        .def("assign", &Surface::operator=)
        //.def(py::self + Surface())
        .def("keep_largest_connected_component", &Surface::keep_largest_connected_component, DOC(Surface, keep_largest_connected_component))
        .def("repair_self_intersections", &Surface::repair_self_intersections, py::arg("volume_threshold") = 0.01, py::arg("cap_threshold") = 170,
             py::arg("needle_threshold") = 2.7, py::arg("collapse_threshold") = 0.14, DOC(Surface, repair_self_intersections))

        .def("remove_small_components", &Surface::remove_small_components, py::arg("volume_threshold") = 30, DOC(Surface, remove_small_components))
        .def("implicit_surface", &Surface::implicit_surface<Surface_implicit_function>, py::arg("implicit_function"), py::arg("bounding_sphere_radius"),
             py::arg("angular_bound") = 30, py::arg("radius_bound") = 0.1, py::arg("distance_bound") = 0.1, DOC(Surface, implicit_surface))

        // TODO : Clean up clip functions
        .def("clip", py::overload_cast<double, double, double, double, bool>(&Surface::clip), py::arg("x0"),
             py::arg("x1"),
             py::arg("x2"),
             py::arg("x3"),
             py::arg("preserve_manifold") = true, DOC(Surface, clip))

        .def("clip", py::overload_cast<Point_3, Vector_3, bool>(&Surface::clip), py::arg("point"), py::arg("vector"), py::arg("preserve_manifold") = true, DOC(Surface, clip, 2))
        .def("clip", py::overload_cast<Plane_3, bool>(&Surface::clip), py::arg("plane"), py::arg("preserve_manifold") = true)
        .def("clip", py::overload_cast<Surface, bool, bool>(&Surface::clip), py::arg("surface"), py::arg("invert") = false, py::arg("preserve_manifold") = true, DOC(Surface, clip, 3))
        .def("clip", py::overload_cast<Point_3, Vector_3, double, bool, bool>(&Surface::clip),
             py::arg("point"), py::arg("vector"), py::arg("radius"), py::arg("invert") = false, py::arg("preserve_manifold") = true, DOC(Surface, clip, 4))

        //.def("clipX", &Surface::clipX ) // TODO: TEST
        .def("get_slice", py::overload_cast<double, double, double, double>(&Surface::get_slice<Slice>), DOC(Surface, get_slice))
        .def("get_slice", py::overload_cast<Plane_3>(&Surface::get_slice<Slice>), DOC(Surface, get_slice, 2))
        //.def("get_slice", py::overload_cast<Point_3,Vector_3>(&Surface::get_slice<Slice>),DOC(Surface,get_slice,3))

        .def("clear", &Surface::clear, DOC(Surface, clear))
        .def("intersection", &Surface::surface_intersection, DOC(Surface, surface_intersection))
        .def("union", &Surface::surface_union, DOC(Surface, surface_union))
        .def("difference", &Surface::surface_difference, DOC(Surface, surface_difference))
        .def("span", &Surface::span, DOC(Surface, span))
        .def("save", &Surface::save, DOC(Surface, save))
        .def("fill_holes", &Surface::fill_holes, DOC(Surface, fill_holes))
        .def("triangulate_faces", &Surface::triangulate_faces, DOC(Surface, triangulate_faces))
        .def("isotropic_remeshing", py::overload_cast<double, unsigned int, bool>(&Surface::isotropic_remeshing), DOC(Surface, isotropic_remeshing))
        .def("adjust_boundary", &Surface::adjust_boundary, DOC(Surface, adjust_boundary))
        .def("smooth_laplacian", &Surface::smooth_laplacian, DOC(Surface, smooth_laplacian))
        .def("smooth_taubin", &Surface::smooth_taubin, DOC(Surface, smooth_taubin))
        .def("smooth_shape", &Surface::smooth_shape, DOC(Surface, smooth_shape))
        .def("make_cube", py::overload_cast<Point_3, Point_3, double>(&Surface::make_cube), DOC(Surface, make_cube))
        .def("make_cube", py::overload_cast<double, double, double, double, double, double, double>(&Surface::make_cube),
             py::arg("x0"), py::arg("y0"), py::arg("z0"), py::arg("x1"), py::arg("y1"), py::arg("z1"), py::arg("edge_length"), DOC(Surface, make_cube, 2))

        .def("make_cone", py::overload_cast<double, double, double, double, double, double, double, double, double>(&Surface::make_cone), DOC(Surface, make_cone))
        .def("make_cone", py::overload_cast<Point_3, Point_3, double, double, double>(&Surface::make_cone), DOC(Surface, make_cone, 2))

        .def("make_cylinder", py::overload_cast<double, double, double, double, double, double, double, double>(&Surface::make_cylinder), DOC(Surface, make_cylinder))
        .def("make_cylinder", py::overload_cast<Point_3, Point_3, double, double>(&Surface::make_cylinder), DOC(Surface, make_cylinder, 2))

        .def("make_sphere", py::overload_cast<double, double, double, double, double>(&Surface::make_sphere), DOC(Surface, make_sphere))
        .def("make_sphere", py::overload_cast<Point_3, double, double>(&Surface::make_sphere), DOC(Surface, make_sphere, 2))

        .def("make_circle_in_plane", py::overload_cast<double, double, double, double, double, double, double, double>(&Surface::make_circle_in_plane))
        .def("make_circle_in_plane", py::overload_cast<Point_3, Vector_3, double, double>(&Surface::make_circle_in_plane))

        .def("is_point_inside", &Surface::is_point_inside, DOC(Surface, is_point_inside))

        .def("get_closest_points", &Surface::get_closest_points, py::arg("point"), py::arg("num") = 8, DOC(Surface, get_closest_points))

        .def("mean_curvature_flow", &Surface::mean_curvature_flow, DOC(Surface, mean_curvature_flow))

        .def("get_shortest_surface_path", py::overload_cast<double, double, double, double, double, double>(&Surface::get_shortest_surface_path), DOC(Surface, get_shortest_surface_path))
        .def("get_shortest_surface_path", py::overload_cast<Point_3, Point_3>(&Surface::get_shortest_surface_path), DOC(Surface, get_shortest_surface_path, 2))

        .def("embed", &Surface::embed, py::arg("other"), py::arg("adjustment") = -0.8, py::arg("smoothing") = 0.4, py::arg("max_iter") = 400, DOC(Surface, embed))
        .def("enclose", &Surface::enclose, py::arg("other"), py::arg("adjustment") = 0.8, py::arg("smoothing") = -0.4, py::arg("max_iter") = 400, DOC(Surface, enclose))
        .def("expose", &Surface::expose, py::arg("other"), py::arg("adjustment") = -0.8, py::arg("smoothing") = 0.4, py::arg("max_iter") = 400, DOC(Surface, expose))
        .def("separate", &Surface::separate, py::arg("other"), py::arg("adjustment") = 0.8, py::arg("smoothing") = 0.4, py::arg("max_iter") = 400, DOC(Surface, separate))

        .def("collapse_edges", py::overload_cast<const double>(&Surface::collapse_edges), DOC(Surface, collapse_edges))
        .def("collapse_edges", py::overload_cast<>(&Surface::collapse_edges), DOC(Surface, collapse_edges, 2))
        .def("split_edges", &Surface::split_edges, DOC(Surface, split_edges))

        // TODO ReMOVE ??
        .def("intersecting_polylines", py::overload_cast<Point_3, Vector_3>(&Surface::polylines_in_plane))
        .def("intersecting_polylines", py::overload_cast<Plane_3>(&Surface::polylines_in_plane))

        .def("extension", py::overload_cast<double, double, double, double, double, double, bool>(&Surface::cylindrical_extension), py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("radius"),
             py::arg("length"),
             py::arg("edge_length"),
             py::arg("use_normal"), DOC(Surface, cylindrical_extension))
        .def("extension", py::overload_cast<const Point_3 &, double, double, double, bool>(&Surface::cylindrical_extension), py::arg("point"),
             py::arg("radius"),
             py::arg("length"),
             py::arg("edge_length"),
             py::arg("use_normal"), DOC(Surface, cylindrical_extension, 2))
        .def("connection", &Surface::cylindrical_connection, DOC(Surface, cylindrical_connection))

        .def("separate_narrow_gaps", &Surface::separate_narrow_gaps, py::arg("adjustment") = -0.5, py::arg("smoothing") = 0.0, py::arg("max_iter") = 400, DOC(Surface, separate_narrow_gaps))
        .def("separate_close_vertices", &Surface::separate_close_vertices, py::arg("adjustment") = 0.5, py::arg("max_iter") = 400, DOC(Surface, separate_close_vertices))

        .def("reconstruct", py::overload_cast<double, double, double>(&Surface::reconstruct),
             py::arg("angular_bound") = 20, py::arg("radius_bound") = 0.1, py::arg("distance_bound") = 0.1, DOC(Surface, reconstruct))

        .def("reconstruct", py::overload_cast<std::string, double, double, double>(&Surface::reconstruct),
             py::arg("filename"), py::arg("angular_bound") = 20, py::arg("radius_bound") = 0.1, py::arg("distance_bound") = 0.1, DOC(Surface, reconstruct, 2))

        .def("convex_hull", &Surface::convex_hull, DOC(Surface, convex_hull))
        .def("num_faces", &Surface::num_faces, DOC(Surface, num_faces))
        .def("num_edges", &Surface::num_edges, DOC(Surface, num_edges))
        .def("num_self_intersections", &Surface::num_self_intersections, DOC(Surface, num_self_intersections))
        .def("num_vertices", &Surface::num_vertices, DOC(Surface, num_vertices))
        .def("distance", &Surface::distance_to_point, DOC(Surface, distance_to_point))
        .def("centeroid", &Surface::centeroid, DOC(Surface, centeroid))
        .def("area", &Surface::area, DOC(Surface, area))
        .def("volume", &Surface::volume, DOC(Surface, volume));

    py::class_<Domain, std::shared_ptr<Domain>>(m, "Domain", DOC(Domain))
        .def(py::init<Surface &, double>(), py::arg("surface"), py::arg("error_bound") = 1.e-7, DOC(Domain, Domain))
        .def(py::init<std::vector<Surface>, double>(), py::arg("surfaces"), py::arg("error_bound") = 1.e-7, DOC(Domain, Domain, 2))
        .def(py::init<std::vector<Surface>, std::shared_ptr<AbstractMap>, double>(), py::arg("surfaces"), py::arg("map"), py::arg("error_bound") = 1.e-7, DOC(Domain, Domain, 3))

        .def("create_mesh", py::overload_cast<double, double, double, double, double, double>(&Domain::create_mesh),
             py::arg("edge_size"), py::arg("cell_size"), py::arg("facet_size"),
             py::arg("facet_angle"), py::arg("facet_distance"), py::arg("cell_radius_edge_ratio"), DOC(Domain, create_mesh))

        .def("create_mesh", py::overload_cast<double>(&Domain::create_mesh), DOC(Domain, create_mesh, 2))
        .def("create_mesh", py::overload_cast<>(&Domain::create_mesh), DOC(Domain, create_mesh, 3))

        .def("radius_ratios_min_max", &Domain::radius_ratios_min_max, DOC(Domain, radius_ratios_min_max))
        .def("dihedral_angles_min_max", &Domain::dihedral_angles_min_max, DOC(Domain, dihedral_angles_min_max))

        .def("radius_ratios", &Domain::radius_ratios, DOC(Domain, radius_ratios))
        .def("dihedral_angles", &Domain::dihedral_angles, DOC(Domain, dihedral_angles))

        .def("get_boundary", &Domain::get_boundary<Surface>, py::arg("tag") = 0, DOC(Domain, get_boundary))
        .def("get_boundaries", &Domain::get_boundaries<Surface>, DOC(Domain, get_boundaries))

        //.def("get_borders", &Domain::get_borders) //TODO DOC(Domain,get_borders)) //TODO
        .def("get_interface", &Domain::get_interface<Surface>) //, DOC(Domain,get_interface))
        .def("get_curve_tags", &Domain::get_curve_tags, DOC(Domain, get_curve_tags))
        .def("get_patches", &Domain::get_patches, DOC(Domain, get_patches))
        .def("get_subdomains", &Domain::get_subdomains, DOC(Domain, get_subdomains))
        .def("check_mesh_connections", &Domain::check_mesh_connections, DOC(Domain, check_mesh_connections))
        .def("lloyd", &Domain::lloyd, py::arg("time_limit") = 0,
             py::arg("max_iter") = 0,
             py::arg("convergence") = 0.02,
             py::arg("freeze_bound") = 0.01,
             py::arg("do_freeze") = true, DOC(Domain, lloyd))

        .def("odt", &Domain::odt, py::arg("time_limit") = 0,
             py::arg("max_iter") = 0,
             py::arg("convergence") = 0.02,
             py::arg("freeze_bound") = 0.01,
             py::arg("do_freeze") = true, DOC(Domain, odt))

        .def("exude", &Domain::exude, py::arg("time_limit") = 0, py::arg("sliver_bound") = 0, DOC(Domain, exude))
        .def("perturb", &Domain::perturb, py::arg("time_limit") = 0, py::arg("sliver_bound") = 0, DOC(Domain, perturb))

        // TODO add sharp border edges multiple surfaces
        .def("add_sharp_border_edges", py::overload_cast<Surface &, double>(&Domain::add_sharp_border_edges<Surface>), py::arg("surface"),
             py::arg("threshold") = 60,
             DOC(Domain, add_sharp_border_edges))
        .def("add_sharp_border_edges", py::overload_cast<Surface &, Plane_3, double>(&Domain::add_sharp_border_edges<Surface, Plane_3>), py::arg("surface"),
             py::arg("plane"),
             py::arg("threshold") = 60,
             DOC(Domain, add_sharp_border_edges, 2))
        .def("clear_borders", &Domain::clear_borders, DOC(Domain, clear_borders))
        .def("clear_features", &Domain::clear_features, DOC(Domain, clear_features))

        .def("remove_subdomain", py::overload_cast<std::vector<int>>(&Domain::remove_subdomain), DOC(Domain, remove_subdomain))
        .def("remove_subdomain", py::overload_cast<int>(&Domain::remove_subdomain), DOC(Domain, remove_subdomain, 2))

        .def("number_of_cells", &Domain::number_of_cells, DOC(Domain, number_of_cells))
        .def("number_of_subdomains", &Domain::number_of_subdomains, DOC(Domain, number_of_subdomains))
        .def("number_of_curves", &Domain::number_of_curves, DOC(Domain, number_of_curves))
        .def("number_of_patches", &Domain::number_of_patches, DOC(Domain, number_of_patches))
        .def("number_of_surfaces", &Domain::number_of_surfaces, DOC(Domain, number_of_surfaces))
        .def("number_of_facets", &Domain::number_of_facets, DOC(Domain, number_of_facets))
        .def("number_of_vertices", &Domain::number_of_vertices, DOC(Domain, number_of_vertices))

        // .def("subdomain_reduction", &Domain::subdomain_reduction<Surface>)

        .def("boundary_segmentations", py::overload_cast<std::pair<int, int>, double>(&Domain::boundary_segmentations<Surface>),
             py::arg("interface"),
             py::arg("angle_in_degree") = 85,
             DOC(Domain, boundary_segmentations))

        .def("boundary_segmentations", py::overload_cast<int, double>(&Domain::boundary_segmentations<Surface>),
             py::arg("subdomain_tag"),
             py::arg("angle_in_degree") = 85,
             DOC(Domain, boundary_segmentations, 2))

        .def("boundary_segmentations", py::overload_cast<double>(&Domain::boundary_segmentations<Surface>),
             py::arg("angle_in_degree") = 85,
             DOC(Domain, boundary_segmentations, 3))

        .def("add_feature", &Domain::add_feature, DOC(Domain, add_feature))
        .def("add_border", &Domain::add_border, DOC(Domain, add_border))
        .def("save", py::overload_cast<std::string, bool>(&Domain::save),
             py::arg("OutPath"),
             py::arg("save_1Dfeatures") = true,
             DOC(Domain, save));

    m.def("load_points", &Wrapper_load_points); // TODO
    m.def("convex_hull", &Wrapper_convex_hull); // TODO

    m.def("separate_overlapping_surfaces", py::overload_cast<Surface &, Surface &, Surface &, double, double, int>(&separate_surface_overlapp<Surface>),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("other"),
          py::arg("edge_movement") = -0.5,
          py::arg("smoothing") = 0.25,
          py::arg("max_iter") = 400,
          DOC(separate_surface_overlapp));

    m.def("separate_overlapping_surfaces", py::overload_cast<Surface &, Surface &, double, double, int>(&separate_surface_overlapp<Surface>),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("edge_movement") = -0.5,
          py::arg("smoothing") = 0.25,
          py::arg("max_iter") = 400,
          DOC(separate_surface_overlapp, 2));

    m.def("separate_close_surfaces", py::overload_cast<Surface &, Surface &, Surface &, double, double, int>(&separate_close_surfaces<Surface>),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("other"),
          py::arg("edge_movement") = -0.5,
          py::arg("smoothing") = 0.25,
          py::arg("max_iter") = 400,
          DOC(separate_close_surfaces));

    m.def("separate_close_surfaces", py::overload_cast<Surface &, Surface &, double, double, int>(&separate_close_surfaces<Surface>),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("edge_movement") = -0.5,
          py::arg("smoothing") = 0.25,
          py::arg("max_iter") = 400,
          DOC(separate_close_surfaces, 2));

    m.def("union_partially_overlapping_surfaces", &union_partially_overlapping_surfaces<Surface>,
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("angle_in_degress") = 36.87,
          py::arg("adjustment") = 0.7,
          py::arg("smoothing") = 0.25,
          py::arg("max_iter") = 8,
          DOC(union_partially_overlapping_surfaces));
}
