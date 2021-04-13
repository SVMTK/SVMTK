#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/complex.h>
#include <pybind11/chrono.h>

#include "Surface.h"
#include "Domain.h"
#include "Slice.h"


namespace py = pybind11;

template <typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

class PyAbstractMap : public AbstractMap{
public:
       using AbstractMap::AbstractMap; /* Inherit constructors */

};

typedef std::function<double(double,double,double)> Surface_implicit_function;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Vector_3 Vector_3;

Vector_3 Wrapper_vector_3(double x,double y,double z)
{
       return Vector_3(x,y,z);
}
Point_3 Wrapper_point_3(double x,double y,double z)
{
       return Point_3(x,y,z);
}
Point_2 Wrapper_point_2(double x,double y)
{
       return Point_2(x,y);
}
Plane_3 Wrapper_plane_3(double x1, double x2 , double x3 , double x4)
{
      return Plane_3(x1,x2,x3,x4);
}
Plane_3 Wrapper_plane_3(Point_3 p1, Vector_3 v1)
{
      return Plane_3(p1,v1);
}
Plane_3 Wrapper_plane_3(Point_3 p1, Point_3 p2 , Point_3 p3)
{
      return Plane_3(p1,p2,p3);
}

std::shared_ptr< Surface > Wrapper_convex_hull(py::array_t< double > point3_array)
{
    py::buffer_info buffer = point3_array.request();

    // Some sanity checks
    if (buffer.ndim != 2)
        throw std::runtime_error("Expected 2d array");

    if (buffer.size % 3 != 0)
        throw std::runtime_error("Point array size must be shape (N, 3)");

    double *ptr_point3_array = (double *) buffer.ptr;

    auto point_vector = std::vector<Surface::Point_3>(buffer.size / 3);

    size_t point_vector_counter = 0;
    for (int i = 0; i < buffer.size; i += 3)
    {
        double x = ptr_point3_array[i];
        double y = ptr_point3_array[i + 1];
        double z = ptr_point3_array[i + 2];

        point_vector[point_vector_counter++] = Point_3(x, y, z);
    }

    auto surface = convex_hull<Surface,Point_3>(point_vector);
    return surface;
}



PYBIND11_MODULE(SVMTK, m) {
   m.doc() = "Surface Volume Meshing Toolkit";
   
   static py::exception<PreconditionError> ex1(m, "PreconditionError");
   static py::exception<EmptyMeshError> ex2(m, "EmptyMeshError");   
   static py::exception<InvalidArgumentError> ex3(m, "InvalidArgumentError");
   static py::exception<AlgorithmError> ex4(m, " AlgorithmError");   
  
   
   py::register_exception_translator([](std::exception_ptr p) 
   {   
       try  {
           if (p) std::rethrow_exception(p);
       } catch (const PreconditionError &e) {
           ex1(e.what());}  
       catch (const InvalidArgumentError &e) {
         ex3(e.what());}   
       catch (const EmptyMeshError &e) {
         ex2(e.what());}   
       catch (const AlgorithmError &e) {
         ex4(e.what());}        
         
   });
   
   
   py::class_<Vector_3,std::shared_ptr<Vector_3>>(m, "Vector_3")
       .def(py::init<double,double,double>())
        .def(py::init<Point_3,Point_3>())
       .def("__repr__",[](Vector_3 const & self)   
        {
           std::ostringstream os;
           os << "(" << self.x() <<", "<< self.y()<<", "<< self.z() <<")";
           return os.str();

        })
       .def(float() * py::self)
       .def(py::self * float())

       
       .def(py::self + py::self)
       .def("squared_length", &Vector_3::squared_length)
       .def("x", &Vector_3::x)
       .def("y", &Vector_3::y)
       .def("z", &Vector_3::z);

   py::class_<Plane_3,std::shared_ptr<Plane_3>>(m, "Plane_3")
       .def(py::init<double,double,double,double>())
       .def(py::init<Point_3,Vector_3>())
       .def("__repr__",[](Plane_3 const & self)   
        {
           std::ostringstream os;
           os << "(" << self.a() <<", "<< self.b()<<", "<< self.c() << "," << self.d() <<")";
           return os.str();

        })
       .def("a", &Plane_3::a)
       .def("b", &Plane_3::b)
       .def("d", &Plane_3::d)
       .def("c", &Plane_3::c);

    py::class_<Point_3,std::shared_ptr<Point_3>>(m, "Point_3") 
       .def(py::init<double,double,double>())
       .def("__repr__",[](Point_3 const & self)   
        {
           std::ostringstream os;
           os << "(" << self.x() <<", "<< self.y()<<", "<< self.z() <<")";
           return os.str();
        })
       .def("__eq__",[](Point_3 const & self, Point_3 const & other)
        {
           return (self.x()==other.x() and self.y()==other.y() and self.z()==other.z()) ; 
        }  ,py::is_operator() ) 
       .def("__ne__",[](Point_3 const & self, Point_3 const & other)
        {
           return (self.x()!=other.x() and self.y()!=other.y() and self.z()!=other.z())   ; 
        }  ,py::is_operator() ) 

       .def("__add__", [](const Point_3 p1, const Vector_3 v1) { p1 + v1; },py::is_operator())
       
       .def(py::self + Vector_3())
       .def(py::self += Vector_3())
       .def(py::self - py::self)
       
       
       
       .def("x", &Point_3::x)
       .def("y", &Point_3::y)
       .def("z", &Point_3::z);
 
    py::class_<Point_2,std::shared_ptr<Point_2>>(m, "Point_2")
       .def(py::init<double,double>())
       .def("__repr__",[](Point_2 const & self)   
        {
           std::ostringstream os;
           os << "(" << self.x() <<", "<< self.y()<<")";
           return os.str();

        })
       .def("x", &Point_2::x)
       .def("y", &Point_2::y);     
 
    py::class_<AbstractMap,PyAbstractMap,std::shared_ptr<AbstractMap>>(m,"AbstractMap");

    py::class_<SubdomainMap,AbstractMap,std::shared_ptr<SubdomainMap>>(m, "SubdomainMap")
        .def(py::init<int>(), py::arg("num_surfaces")=0)
        .def("print",  &SubdomainMap::print)
        .def("set_number_of_surfaces", &SubdomainMap::set_number_of_surfaces)
        .def("add_interface", &SubdomainMap::add_interface)
        .def("get_interfaces", &SubdomainMap::get_interfaces)        
        .def("erase", &SubdomainMap::erase)
        .def("get_map", &SubdomainMap::get_map)
        .def("get_tags", &SubdomainMap::get_tags) 
        .def("add", &SubdomainMap::add);
       

    py::class_<Slice,std::shared_ptr<Slice>>(m, "Slice")
        .def(py::init<>())
        .def(py::init<Slice&>())
        .def(py::init<Plane_3>())
        .def(py::init<double,double,double,double>()) 
        .def("create_mesh", &Slice::create_mesh) 
        .def("simplify", &Slice::simplify) 
        .def("save", &Slice::save)
        .def("slice_surfaces", &Slice::slice_surfaces<Surface> ) 
        .def("as_surface", &Slice::as_surface<Surface>  ) 
        .def("add_surface_domains", py::overload_cast<std::vector<Surface>, AbstractMap&>( &Slice::add_surface_domains<Surface> ) ) 
        .def("add_surface_domains", py::overload_cast<std::vector<Surface>>( &Slice::add_surface_domains<Surface> ) ) 
        .def("number_of_constraints",&Slice::number_of_constraints)
        .def("number_of_subdomains",&Slice::number_of_subdomains)
        .def("number_of_faces",&Slice::number_of_faces)
        .def("connected_components",&Slice::connected_components) 
        .def("keep_largest_connected_component",&Slice::keep_largest_connected_component)
        .def("remove_subdomain", py::overload_cast<int> ( &Slice::remove_subdomain)) 
        .def("remove_subdomain", py::overload_cast<std::vector<int>> ( &Slice::remove_subdomain)) 
        .def("get_constraints", &Slice::get_constraints )
        .def("add_constraint", &Slice::add_constraint )
        .def("add_constraints",py::overload_cast<Slice&>( &Slice::add_constraints));


    py::class_<Surface,std::shared_ptr<Surface>>(m, "Surface")
        .def(py::init<std::string &>())
        .def(py::init<>())
        .def(py::init<Surface&>())     
        
        
        .def("__copy__",  [](const Surface &self) 
         {return Surface(self);})
        .def("copy",  [](const Surface &self) 
         {return Surface(self);})
        .def("assign", &Surface::operator=)
        .def("keep_largest_connected_component",&Surface::keep_largest_connected_component)
        
        .def("implicit_surface",  &Surface::implicit_surface<Surface_implicit_function>, py::arg("implicit_function") ,py::arg("bounding_sphere_radius"),
                                                             py::arg("angular_bound")=30,py::arg("radius_bound")=0.1,py::arg("distance_bound")=0.1)

        .def("clip", py::overload_cast<double,double,double,double,bool>( &Surface::clip ), py::arg("x0"),py::arg("x1"),py::arg("x2"),py::arg("x3"),py::arg("preserve_manifold")=true ) 
        .def("clip",py::overload_cast<Point_3,Vector_3,bool>( &Surface::clip ), py::arg("point"),py::arg("vector"),py::arg("preserve_manifold")=true )         
        .def("clip",py::overload_cast<Plane_3,bool>( &Surface::clip ), py::arg("plane"),py::arg("preserve_manifold")=true )        
        .def("clip",py::overload_cast<Surface,bool,bool>( &Surface::clip ), py::arg("surface"), py::arg("invert")=false,py::arg("preserve_manifold")=true )  
        .def("clip",py::overload_cast<Point_3,Vector_3,double,bool,bool>( &Surface::clip ),
                                             py::arg("point"),py::arg("vector"),py::arg("radius"),py::arg("invert")=false,py::arg("preserve_manifold")=true )            
        
        .def("slice", py::overload_cast<double , double, double , double>(&Surface::mesh_slice<Slice>)) 

        .def("clear" , &Surface::clear) 
        .def("intersection", &Surface::surface_intersection)
        .def("union", &Surface::surface_union)
        .def("difference", &Surface::surface_difference)

        .def("span", &Surface::span) 
        .def("save", &Surface::save)

        .def("fill_holes", &Surface::fill_holes)
        .def("triangulate_faces", &Surface::triangulate_faces)
        .def("isotropic_remeshing", py::overload_cast<double , unsigned int , bool>(&Surface::isotropic_remeshing))
        .def("adjust_boundary", &Surface::adjust_boundary)

        .def("smooth_laplacian", &Surface::smooth_laplacian)
        .def("smooth_taubin", &Surface::smooth_taubin)
        .def("smooth_shape", &Surface::smooth_shape)

        .def("make_cube", py::overload_cast<Point_3,Point_3,double>( &Surface::make_cube))
        .def("make_cube", py::overload_cast<double,double,double,double,double,double,double>(&Surface::make_cube),
                          py::arg("x0"),py::arg("y0"),py::arg("z0"),py::arg("x1"),py::arg("y1"),py::arg("z1"),  py::arg("N")=10)
                          
	.def("make_cone", py::overload_cast<double,double,double,double,double,double,double,double,double>(&Surface::make_cone))
	.def("make_cone", py::overload_cast<Point_3,Point_3,double,double,double>(&Surface::make_cone))
	.def("make_cylinder", py::overload_cast<double,double,double,double,double,double,double,double>(&Surface::make_cylinder))	
	.def("make_cylinder", py::overload_cast<Point_3,Point_3,double,double>(&Surface::make_cylinder))
        .def("make_sphere", py::overload_cast<double,double,double,double,double>(&Surface::make_sphere ))
        .def("make_sphere", py::overload_cast<Point_3,double,double>(&Surface::make_sphere ))
        .def("make_circle_in_plane",  py::overload_cast<double,double,double,double,double,double,double,double>(&Surface::make_circle_in_plane))
        .def("make_circle_in_plane",  py::overload_cast<Point_3,Vector_3,double,double>(&Surface::make_circle_in_plane))
        
        .def("does_bound_volume", &Surface::does_bound_a_volume)
        .def("is_point_inside", &Surface::is_point_inside)
        .def("get_closest_points", &Surface::get_closest_points, py::arg("p1"),py::arg("num")=8)

        .def("mean_curvature_flow", &Surface::mean_curvature_flow)
        .def("get_shortest_surface_path", py::overload_cast<double , double, double , double,double,double>( &Surface::get_shortest_surface_path) )
        .def("get_shortest_surface_path", py::overload_cast<Point_3,Point_3>( &Surface::get_shortest_surface_path) )

        .def("embed", &Surface::embed, py::arg("other") , py::arg("adjustment")=-0.5,py::arg("smoothing")=0.4, py::arg("max_iter")=400)
        .def("enclose", &Surface::enclose,py::arg("other") , py::arg("adjustment")=0.,py::arg("smoothing")=0.4, py::arg("max_iter")=400)
        .def("separate", &Surface::separate, py::arg("other") , py::arg("adjustment")=-0.5,py::arg("smoothing")=0.4, py::arg("max_iter")=400)


        .def("collapse_edges", py::overload_cast<const double >( &Surface::collapse_edges))
        .def("collapse_edges", py::overload_cast<>(&Surface::collapse_edges))
        .def("split_edges", &Surface::split_edges)
        .def("assert_non_empty_mesh", &Surface::assert_non_empty_mesh)

        .def("extension", py::overload_cast<double,double,double,double,double,double,bool>(&Surface::cylindric_extension) )
        .def("extension", py::overload_cast<const Point_3&, double,double,double,bool>( &Surface::cylindric_extension) )
        .def("connection", &Surface::cylindric_connection)
        .def("separate_narrow_gaps", &Surface::separate_narrow_gaps, py::arg("adjustment")=-0.5,py::arg("smoothing")=0.4,py::arg("max_iter")=400) //TODO test 
        .def("reconstruct", &Surface::reconstruct, py::arg("angular_bound")=20,py::arg("radius_bound")=0.1,py::arg("distance_bound")=0.1 )
        .def("convex_hull", &Surface::convex_hull)

        .def("num_faces", &Surface::num_faces)
        .def("num_edges", &Surface::num_edges)
        .def("num_self_intersections", &Surface::num_self_intersections) 
        .def("num_vertices", &Surface::num_vertices)
        .def("distance", &Surface::distance_to_point) 
        .def("centeroid", &Surface::centeroid)
        .def("area", &Surface::area)
        .def("volume", &Surface::volume);    


    py::class_<Domain,std::shared_ptr<Domain>>(m, "Domain")
        .def(py::init<Surface &,double>(), py::arg("surface"), py::arg("accuracy")=1.e-7)
        .def(py::init<std::vector<Surface>,double>(),py::arg("surfaces"), py::arg("accuracy")=1.e-7)
        .def(py::init<std::vector<Surface>, std::shared_ptr<AbstractMap>,double>(), py::arg("surfaces"), py::arg("map"), py::arg("accuracy")=1.e-7)

        .def("create_mesh", py::overload_cast<double,double,double,double,double,double>( &Domain::create_mesh), 
                            py::arg("edge_size"), py::arg("cell_size"), py::arg("facet_size"),
                            py::arg("facet_angle"),py::arg("facet_distance"), py::arg("cell_radius_edge_ratio") )  

        .def("create_mesh", py::overload_cast<double>(&Domain::create_mesh))
        .def("radius_ratio_min_max", &Domain::radius_ratio_min_max)
        .def("dihedral_angles_min_max", &Domain::dihedral_angles_min_max)
        .def("radius_ratio", &Domain::radius_ratio)
        .def("dihedral_angles", &Domain::dihedral_angles)       

        .def("get_boundary", &Domain::get_boundary<Surface>, py::arg("tag")=0)
        .def("get_boundaries", &Domain::get_boundaries<Surface>)
         .def("get_borders", &Domain::get_borders)
        .def("get_curves", &Domain::get_curves)
        .def("get_patches", &Domain::get_patches)
        .def("get_subdomains", &Domain::get_subdomains)

        .def("lloyd", &Domain::lloyd,     py::arg("time_limit")=0, py::arg("max_iteration_number")=0, py::arg("convergence")=0.02, py::arg("freeze_bound")=0.01,py::arg("do_freeze")=true)
        .def("odt", &Domain::odt,         py::arg("time_limit")=0, py::arg("max_iteration_number")=0, py::arg("convergence")=0.02, py::arg("freeze_bound")=0.01,py::arg("do_freeze")=true)
        .def("exude", &Domain::exude,     py::arg("time_limit")=0, py::arg("sliver_bound")=0)
        .def("perturb", &Domain::perturb, py::arg("time_limit")=0, py::arg("sliver_bound")=0)

        .def("add_sharp_border_edges", py::overload_cast<Surface&,double>( &Domain::add_sharp_border_edges<Surface>), py::arg("surface") , py::arg("threshold")=60 )
        .def("clear_borders", &Domain::clear_borders)
        .def("clear_features", &Domain::clear_features)
        .def("remove_subdomain", py::overload_cast<std::vector<int>>(&Domain::remove_subdomain))
        .def("remove_subdomain", py::overload_cast<int>(&Domain::remove_subdomain))
 
        .def("number_of_cells", &Domain::number_of_cells)
        .def("number_of_subdomains", &Domain::number_of_subdomains) 
        .def("number_of_curves", &Domain::number_of_curves)
        .def("number_of_patches", &Domain::number_of_patches)
        .def("number_of_surfaces", &Domain::number_of_surfaces)
        .def("number_of_facets", &Domain::number_of_facets)
        .def("number_of_vertices", &Domain::number_of_vertices)

        .def("boundary_segmentations",py::overload_cast<int,double>( &Domain::boundary_segmentations<Surface>),py::arg("subdomain_tag") , py::arg("angle_in_degree")=85  )
        .def("boundary_segmentations",py::overload_cast<double>(&Domain::boundary_segmentations<Surface> ),py::arg("angle_in_degree")=85 )
        .def("add_feature", &Domain::add_feature) 
        .def("add_border", &Domain::add_border)
        .def("save", &Domain::save, py::arg("OutPath"), py::arg("save_1Dfeatures")=true); 
                  


       m.def("convex_hull", &Wrapper_convex_hull); 
       //TODO : Rename edge_movement.
       m.def("separate_overlapping_surfaces",  py::overload_cast<Surface&,Surface&,Surface&,double,double,int>( &separate_surface_overlapp<Surface>),
                                   py::arg("surf1"), py::arg("surf2"), py::arg("other"), 
                                   py::arg("edge_movement")=-0.25 , py::arg("smoothing")=0.3 ,py::arg("max_iter")=400  );

       m.def("separate_overlapping_surfaces",  py::overload_cast<Surface&,Surface&,double,double,int>( &separate_surface_overlapp<Surface>) ,
                                   py::arg("surf1"), py::arg("surf2"), 
                                   py::arg("edge_movement")=-0.25 , py::arg("smoothing")=0.3, py::arg("max_iter")=400  );

       m.def("separate_close_surfaces",  py::overload_cast<Surface&,Surface&,Surface&,double,double,int>( &separate_close_surfaces<Surface>) ,
                                   py::arg("surf1"), py::arg("surf2"), py::arg("other"), 
                                   py::arg("edge_movement")=-0.25 , py::arg("smoothing")=0.3,py::arg("max_iter")=400 );

       m.def("separate_close_surfaces",  py::overload_cast<Surface&,Surface&,double,double,int>(&separate_close_surfaces<Surface>) ,
                                   py::arg("surf1"), py::arg("surf2"), 
                                   py::arg("edge_movement")=-0.25 , py::arg("smoothing")=0.3 ,py::arg("max_iter")=400 );

       m.def("union_partially_overlapping_surfaces", &union_partially_overlapping_surfaces<Surface>,
                                   py::arg("surf1"), py::arg("surf2"), py::arg("clusterth")=36.87 ,
                                   py::arg("edge_movement")=0.25 , py::arg("smoothing")=1, py::arg("max_iter")=8  ); 



}
