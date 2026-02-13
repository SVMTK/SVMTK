#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>

#include "docstrings.h"
#include "Surface.h"
#include "Domain.h"
#include "Slice.h"

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
typedef Kernel::Vector_2 Vector_2;

/* -- CGAL wrapping */

Vector_3 Wrapper_vector_3(double x, double y, double z)
{
     return Vector_3(x, y, z);
}
Vector_3 Wrapper_vector_3(std::tuple<double,double,double> vector )
{
     return Vector_3(std::get<0>(vector), std::get<1>(vector), std::get<2>(vector));
}

Point_3 Wrapper_point_3(double x, double y, double z)
{
     return Point_3(x, y, z);
}
Point_3 Wrapper_point_3(std::tuple<double,double,double> point )
{
     return Point_3(std::get<0>(point), std::get<1>(point), std::get<2>(point));
}

Vector_2 Wrapper_vector_2(double x, double y)
{
     return Vector_2(x, y);
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
Plane_3 Wrapper_plane_3(std::tuple<double,double, double,double> plane)
{
     return Plane_3(std::get<0>(plane), std::get<1>(plane), std::get<2>(plane), std::get<3>(plane));
}

std::vector<std::array<double,3>> wrapper_point_vector( std::vector<Point_3> cgalpoints )
{    
     std::vector< std::array<double, 3> > points; 
     std::array<double, 3> point;
     for (  auto p : cgalpoints ) 
     {
          point[0] = p.x();
          point[1] = p.y();
          point[2] = p.z();
          points.push_back(point);
     }
     return points;
}

std::vector<std::array<double,2>> wrapper_point_vector( std::vector<Point_2> cgalpoints )
{    
     std::vector< std::array<double, 2> > points; 
     std::array<double, 2> point;
     for (  auto p : cgalpoints ) 
     {
          point[0] = p.x();
          point[1] = p.y();
          points.push_back(point);
     }
     return points;
}


/**
* @brief A collection of overloaded functions that saves object directly into a npz file.
*         Surface  : points, faces, face_markers 
*         Domain   : points, cells, faces ,edges. cell_markers, face_markers, edge_markers
*         Slice    : points, cells, faces, cell_markers, face_markers,
* TODO: remove templates 
*/
template< typename Point, typename Face=std::array<std::size_t,3> >
void numpy_savez_wrapper(std::string& filename, std::vector<Point> points, std::vector<Face> faces, std::vector<int> face_markers) 
{
     py::module_ numpy_module = py::module_::import("numpy");
     py::object save_func = numpy_module.attr("savez");     
     save_func( filename, py::arg("points")=points, py::arg("faces")=faces, py::arg("face_markers")=face_markers );
}

template< typename Point, typename Cell>
void numpy_savez_wrapper(const std::string& filename, std::vector<Point> points, std::vector<Cell> cells)
{ 
     py::module_ numpy_module = py::module_::import("numpy");
     py::object save_func = numpy_module.attr("savez");     
     save_func(filename, py::arg("points")=points, py::arg("cells")=cells);
}

template< typename Point, typename Face, typename Cell>
void numpy_savez_wrapper(const std::string& filename, std::vector<Point> points, 
     std::vector<Cell> cells, std::vector<int> cell_markers ,
     std::vector<Face> faces, std::vector<int> face_markers)
{ 
     py::module_ numpy_module = py::module_::import("numpy");
     py::object save_func = numpy_module.attr("savez");     
     save_func(filename, py::arg("points")=points, 
          py::arg("cells")=cells,py::arg("cell_markers")=cell_markers, 
          py::arg("faces")=faces, py::arg("face_markers")=face_markers );
}

template< typename Point, typename Face, typename Cell, typename Edge>
void numpy_savez_wrapper(const std::string& filename, std::vector<Point> points, 
     std::vector<Cell> cells, std::vector<int> cell_markers ,
     std::vector<Face> faces, std::vector<int> face_markers,
     std::vector<Edge> edges, std::vector<int> edge_markers )
{ 
     py::module_ numpy_module = py::module_::import("numpy");
     py::object save_func = numpy_module.attr("savez");     
     save_func(filename, py::arg("points")=points, 
          py::arg("cells")=cells,py::arg("cell_markers")=cell_markers, 
          py::arg("faces")=faces, py::arg("face_markers")=face_markers,
          py::arg("edges")=faces, py::arg("edge_markers")=edge_markers );
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

     auto surface = convex_hull(point_vector);
     return surface;
}

PYBIND11_MODULE(SVMTK, m)
{
     m.doc() = "Surface Volume Meshing Toolkit";

     PYBIND11_CONSTINIT static py::gil_safe_call_once_and_store<py::object> exc_1_storage;
     exc_1_storage.call_once_and_store_result(
     [&]() { return py::exception<PreconditionError>(m, "PreconditionError"); });
    
     PYBIND11_CONSTINIT static py::gil_safe_call_once_and_store<py::object> exc_2_storage;
     exc_2_storage.call_once_and_store_result(
     [&]() { return py::exception<EmptyMeshError>(m, "EmptyMeshError");  });    
    
     PYBIND11_CONSTINIT static py::gil_safe_call_once_and_store<py::object> exc_3_storage;
     exc_3_storage.call_once_and_store_result(
     [&]() { return py::exception<InvalidArgumentError>(m, "InvalidArgumentError");  });    
    
     py::register_exception_translator([](std::exception_ptr p)
     {
          try
          {
               if( p )
                    std::rethrow_exception(p);
          }    
          catch (const PreconditionError &e)
          {
               py::set_error(exc_1_storage.get_stored(), e.what());
          }
          catch (const EmptyMeshError &e)
          {
               py::set_error(exc_2_storage.get_stored(), e.what());
          }
          catch (const InvalidArgumentError &e)
          {
               py::set_error(exc_3_storage.get_stored(), e.what());
          } 
                              
     });
 
     py::class_<Vector_3, std::shared_ptr<Vector_3>>(m, "Vector_3", DOC(Vector3))
          .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"), DOC(Vector3, Vector3))
          .def(py::init<Point_3, Point_3>(), py::arg("source"), py::arg("target"), DOC(Vector3, Vector3, 2))
          .def(py::init([](std::tuple<double,double,double> vector ){ 
                return Vector_3(  std::get<0>(vector) , std::get<1>(vector), std::get<2>(vector));
           }))

          .def("__repr__", [](Vector_3 const &self)
          {
             std::ostringstream os;
             os << "(" << self.x() << ", " << self.y() << ", " << self.z() << ")";
             return os.str(); 
          })
                 
          .def(float() * py::self)
          .def(py::self * float())
          .def(py::self + py::self)
          .def("squared_length", &Vector_3::squared_length, "Returns the squared length of the Vector")
        
          .def("x", &Vector_3::x, "Returns x direction.")
          .def("y", &Vector_3::y, "Returns y direction.")
          .def("z", &Vector_3::z, "Returns z direction.");

     py::class_<Vector_2, std::shared_ptr<Vector_2>>(m, "Vector_2", DOC(Vector2)) 
          .def(py::init<double, double>(), py::arg("x"), py::arg("y"), DOC(Vector2, Vector2))
          .def(py::init<Point_2, Point_2>(), py::arg("source"), py::arg("target"), DOC(Vector2, Vector2, 2)) 
          .def( py::init( [] (std::tuple<double,double> point ){   
               return Vector_2( std::get<0>( point), std::get<1>(point));
          }))

          .def("__repr__", [](Vector_2 const &self)
          {
               std::ostringstream os;
               os << "(" << self.x() << ", " << self.y() << ")";
               return os.str(); 
          })
             
          .def(float() * py::self)
          .def(py::self * float())
          .def(py::self + py::self)
          .def("squared_length", &Vector_2::squared_length, "Returns the squared length of the Vector")
          .def("x", &Vector_2::x, "Returns x direction.")
          .def("y", &Vector_2::y, "Returns y direction.");
        
     py::class_<Plane_3, std::shared_ptr<Plane_3>>(m, "Plane_3", DOC(Plane3))
          .def(py::init<double, double, double, double>(), py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"), DOC(Plane3, Plane3))
          .def(py::init<Point_3, Vector_3>(), py::arg("point"), py::arg("vector"), DOC(Plane3, Plane3, 2))
          .def(py::init<Point_3, Point_3, Point_3>(), py::arg("p1"), py::arg("p2"), py::arg("p3")) 
          .def( py::init([] (std::tuple<double,double,double, double> plane ){ 
               return Plane_3(  std::get<0>(plane), std::get<1>(plane), std::get<2>(plane), std::get<3>(plane));
          }))

          .def("__repr__", [](Plane_3 const &self)
          {
               std::ostringstream os;
               os << "(" << self.a() << ", " << self.b() << ", " << self.c() << "," << self.d() << ")";
               return os.str(); 
          })
             
          .def("a", &Plane_3::a, "Returns first coefficient of plane equation.")
          .def("b", &Plane_3::b, "Returns second coefficient of plane equation.")
          .def("c", &Plane_3::c, "Returns third coefficient of plane equation.")
          .def("d", &Plane_3::d, "Returns fourth coefficient of plane equation.");

     py::class_<Point_3, std::shared_ptr<Point_3>>(m, "Point_3", DOC(Point3))
          .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"), DOC(Point3, Point3))
          .def( py::init([](std::tuple<double,double,double> point){ 
               return Point_3( std::get<0>(point), std::get<1>(point), std::get<2>(point) );
          }))
          .def("__repr__", [](Point_3 const &self){
             std::ostringstream os;
             os << "(" << self.x() <<", "<< self.y()<<", "<< self.z() <<")";
             return os.str(); 
          })
             
          .def("__attr__",  [](Point_3 const &self) 
               { return ( self.x() , self.y() , self.z() ); })   
           
          .def("__eq__", [](Point_3 const &self, Point_3 const &other)
               { return (self.x() == other.x() and self.y() == other.y() and self.z() == other.z()); }, py::is_operator())
        
          .def("__ne__", [](Point_3 const &self, Point_3 const &other)
            { return (self.x() != other.x() and self.y() != other.y() and self.z() != other.z()); }, py::is_operator())
            
          .def(py::hash(py::self))    
          .def(py::self + Vector_3())
          .def(py::self += Vector_3())
          .def(py::self -= Vector_3()) 
          .def(py::self - Vector_3())
          .def("x", &Point_3::x, "Returns x coordinate.")
          .def("y", &Point_3::y, "Returns y coordinate.")
          .def("z", &Point_3::z, "Returns z coordinate.");

     py::class_<Point_2, std::shared_ptr<Point_2>>(m, "Point_2", DOC(Point2))
          .def(py::init<double, double>(), DOC(Point2, Point2))
          .def(py::init( [] (std::tuple<double,double> point){ 
               return Point_2( std::get<0>( point), std::get<1>(point));
          }))
           
          .def("__attr__",  [](Point_2 const &self) 
               { return ( self.x() , self.y() ); })   
            
          .def(py::hash(py::self))
   
          .def("__eq__", [](Point_2 const &self, Point_2 const &other)
               {return (self.x() == other.x() and  self.y() == other.y()); }, py::is_operator())
            
          .def("__ne__", [](Point_2 const &self, Point_2 const &other)
               {return (self.x() != other.x() and self.y() != other.y()); }, py::is_operator())      
                               
          .def("__repr__", [](Point_2 const &self)
            {
              std::ostringstream os;
              os << "(" << self.x() << ", " << self.y() << ")";
              return os.str(); 
          })

          .def(py::self + Vector_2())
          .def(py::self += Vector_2())
          .def(py::self - Vector_2())                 
          .def(py::self - Point_2())       
          .def("x", &Point_2::x, "Returns x coordinate.")
          .def("y", &Point_2::y, "Returns y coordinate.");

     py::class_<AbstractMap, PyAbstractMap, std::shared_ptr<AbstractMap>>(m, "AbstractMap");

     py::class_<SubdomainMap, AbstractMap, std::shared_ptr<SubdomainMap>>(m, "SubdomainMap", DOC(SubdomainMap))
          .def(py::init<int>(), py::arg("num_surfaces") = 0, DOC(SubdomainMap, SubdomainMap))
   
          .def("__repr__", [](SubdomainMap &self)
               { return self.print(); })

          .def("add",              &SubdomainMap::add,                DOC(SubdomainMap, add))
          .def("print",            &SubdomainMap::print,              DOC(SubdomainMap, print))
          .def("erase",            &SubdomainMap::erase,              DOC(SubdomainMap, erase))
          .def("get_map",          &SubdomainMap::get_map,            DOC(SubdomainMap, get_map))
          .def("get_tags",         &SubdomainMap::get_tags,           DOC(SubdomainMap, get_tags))
          .def("add_interface",    &SubdomainMap::add_interface,      DOC(SubdomainMap, add_interface))
          .def("get_interfaces",   &SubdomainMap::get_interfaces,     DOC(SubdomainMap, get_interfaces))
          .def("make_interfaces",  &SubdomainMap::make_interfaces,    DOC(SubdomainMap, make_interfaces))
          .def("set_num_surfaces", &SubdomainMap::set_num_surfaces,   DOC(SubdomainMap, set_num_surfaces));

     py::class_<Slice, std::shared_ptr<Slice>>(m, "Slice", DOC(Slice))
          .def(py::init<>(), DOC(Slice, Slice))
          .def(py::init<Plane_3>(), DOC(Slice, Slice, 2))
          .def(py::init<Point_3, Vector_3>(), DOC(Slice, Slice, 3))
          .def(py::init<double, double, double, double>(), DOC(Slice, Slice, 4))
          .def(py::init( [] ( std::tuple<double,double,double,double>  plane ) {
               return Slice(Wrapper_plane_3( plane ));
          }), DOC(Slice,Slice,5))

          .def("get_points",[]( Slice &self){
               auto points = self.get_points();
               return py::array(py::cast(std::move(points))); },  py::return_value_policy::copy, DOC(Slice,get_points) )
                            
          .def("get_cells", []( Slice &self){
               auto cells  = self.get_cells();
               return py::array(py::cast(std::move(cells))); }, py::return_value_policy::move ,DOC(Slice,get_cells))  

          .def("get_cell_tags", []( Slice &self){
               auto cell_tags  = self.get_cell_tags();
               return py::array(py::cast(std::move(cell_tags))); }, py::return_value_policy::move , DOC(Slice,get_cell_tags) )          
        
          .def("get_facets", []( Slice &self, bool exclude_unmarked){
               auto facets  = self.get_facets(exclude_unmarked);
               return py::array(py::cast(std::move(facets))); },
               py::arg("exclude_unmarked")=true  ,py::return_value_policy::move , DOC(Slice,get_facets))          
                           
          .def("get_facet_tags", []( Slice &self, bool exclude_unmarked){
               auto facet_tags = self.get_facet_tags(exclude_unmarked);
               return py::array(py::cast(std::move(facet_tags))); },
               py::arg("exclude_unmarked")=true,
               py::return_value_policy::move, DOC(Slice,get_facet_tags) )          

          .def("add_polygon_domain", py::overload_cast<std::vector<Point_2>, int, int >(&Slice::add_polygon_domain) , 
               py::arg("polygon"), 
               py::arg("polygon_tag"), 
               py::arg("replace_tag") = 0, DOC(Slice,add_polygon_domain) )    
 
          
          .def("add_polygon_domain", [](Slice &self, std::vector<std::tuple<double, double>> pyarray, int polygon_tag, int replace_tag ) {
               std::vector<Point_2> wrap; 
               for( auto tpl : pyarray )
                    wrap.emplace_back(  std::get<0>(tpl), std::get<1>(tpl) );
               self.add_polygon_domain( wrap, polygon_tag ,  replace_tag);
          },  DOC(Slice,add_polygon_domain,2)  )

          .def("get_constraints",       &Slice::get_constraints,      DOC(Slice, get_constraints))
          .def("get_feature_vertices",  &Slice::get_feature_vertices, py::return_value_policy::move, DOC(Slice,get_feature_vertices) ) // Check return policy 
              
          .def("create_mesh", py::overload_cast<double>(&Slice::create_mesh),             DOC(Slice, create_mesh))
          .def("create_mesh", py::overload_cast<double, double>(&Slice::create_mesh),     DOC(Slice, create_mesh, 2))
        
          .def("remove_constraints",  py::overload_cast<std::vector<std::vector<Point_2>>>(&Slice::remove_constraints),  DOC(Slice,remove_constraints))      
          .def("remove_constraints",  py::overload_cast<std::vector<Point_2>>(&Slice::remove_constraints),DOC(Slice,remove_constraints))
       
          .def("remove_constraints",  [](Slice &self, std::vector<std::array<double,2>> pyarray) {
               std::vector<Point_2> wrap; 
               for( auto point : pyarray )
                    wrap.emplace_back(std::get<0>(point), std::get<1>(point) );
               self.remove_constraints(wrap);
          }, DOC(Slice,remove_constraints,4))     

          .def("remove_constraints",  [](Slice &self, std::vector<std::vector<std::array<double,2>>> pyarray) {
               for( auto polygon : pyarray)
               {
                    std::vector<Point_2> wrap; 
                    for( auto point : polygon )
                         wrap.emplace_back(std::get<0>(point), std::get<1>(point) );
                    self.remove_constraints(wrap);
               }
          }, DOC(Slice,remove_constraints,3))       

          .def("simplify",    &Slice::simplify,   DOC(Slice, simplify))
        
          .def("num_constraints",  &Slice::num_constraints, DOC(Slice, num_constraints))
          .def("num_subdomains",   &Slice::num_subdomains,  DOC(Slice, num_subdomains))
          .def("num_cells",        &Slice::num_cells,       DOC(Slice, num_cells))
          .def("num_facets",       &Slice::num_facets,      DOC(Slice, num_facets))
          .def("num_vertices",     &Slice::num_vertices,    DOC(Slice, num_vertices))

          .def("slice_mesh", &Slice::slice_mesh<Domain,Surface> , DOC(Slice, slice_mesh))
          .def("slice_surfaces", &Slice::slice_surfaces<Surface>, DOC(Slice, slice_surfaces))
       
          .def("smooth_constraints",    &Slice::smooth_constraints ,  DOC(Slice, smooth_constraints))
          .def("bifurcation_split",     &Slice::bifurcation_split,    DOC(Slice,bifurcation_split) ) 
          .def("remove_small_polygons", &Slice::remove_small_polygons,DOC(Slice,remove_polygons))

          .def("connected_components",             &Slice::connected_components,            DOC(Slice, connected_components))
          .def("keep_largest_connected_component", &Slice::keep_largest_connected_component,DOC(Slice, keep_largest_connected_component))

          .def("remove_subdomain", py::overload_cast<int>(&Slice::remove_subdomain),             DOC(Slice, remove_subdomain))
          .def("remove_subdomain", py::overload_cast<std::vector<int>>(&Slice::remove_subdomain),DOC(Slice, remove_subdomain, 2))

          .def("export_as_surface", py::overload_cast<>(&Slice::export_as_surface<Surface>),                    DOC(Slice, export_as_surface))
          .def("export_as_surface", py::overload_cast<std::vector<double>>(&Slice::export_as_surface<Surface>), DOC(Slice, export_as_surface))  

          .def("reconstruct_constraints",&Slice::reconstruct_constraints, py::arg("target_edge_length"), py::arg("merge")=true, DOC(Slice,reconstruct_constraints))

          .def("add_surface_domains", py::overload_cast<std::vector<Surface>, AbstractMap &>(&Slice::add_surface_domains<Surface>),    DOC(Slice, add_surface_domains))
          .def("add_surface_domains", py::overload_cast<std::vector<Surface>>(&Slice::add_surface_domains<Surface>),                   DOC(Slice, add_surface_domains, 2))
        
          .def("add_constraint", py::overload_cast<std::vector<Point_2>>(&Slice::add_constraint), DOC(Slice, add_constraint))
          .def("add_constraint", [](Slice &self, std::vector<std::array<double,2>> pyarray) {
               std::vector<Point_2> wrap; 
               for( auto tpl : pyarray )
                    wrap.emplace_back(std::get<0>(tpl), std::get<1>(tpl));
               self.add_constraint( wrap);
          } , DOC(Slice, add_constraint,2) ) 

          .def("add_constraints", py::overload_cast<Slice &>(&Slice::add_constraints), DOC(Slice, get_constraints))          
          .def("add_constraints", [](Slice &self, std::vector<std::vector<std::array<double,2>>> pyarray) {
               for( auto constraint : pyarray )
               {
                    std::vector<Point_2> wrap; 
                    for( auto point: constraint )
                         wrap.emplace_back(std::get<0>(point), std::get<1>(point));
                    self.add_constraint(wrap);
               }
          }, DOC(Slice, add_constraints,3))
          
          .def("save", []( Slice &self, std::string filename, bool exclude_unmarked){
               
               self.assert_non_empty_mesh();  
               std::string extension = filename.substr(filename.find_last_of(".")+1);
               std::ofstream out(filename);
               if( extension =="npz")
               { 
                    auto points     = self.get_points();
                    auto cells      = self.get_cells();
                    auto cell_tags  = self.get_cell_tags();
                    auto facets     = self.get_facets(exclude_unmarked);
                    auto facet_tags = self.get_facet_tags(exclude_unmarked);
                    
                    numpy_savez_wrapper(filename, points, cells, cell_tags, facets, facet_tags);
               }
               else
                    self.save(filename);
        
          }, py::arg("filename"), py::arg("exclude_unmarked")=true ,DOC(Slice, save));

 
     py::class_<Surface, std::shared_ptr<Surface>>(m, "Surface", DOC(Surface))
          .def(py::init<>(), DOC(Surface, Surface))
        
          .def(py::init<std::string &>(), py::arg("filename"), DOC(Surface, Surface, 2))
          .def(py::init<Surface &>(),     py::arg("surf"),     DOC(Surface, Surface, 3))     

          .def(py::init<std::vector<Point_3>, int> (), 
               py::arg("points"),  
               py::arg("num_neighbors"), DOC(Surface,Surface,4 ))

          .def( py::init( []( std::vector< std::array<double, 3> > points, int num_neighbors ){
               std::vector<Point_3> cgalpoints;
               for( auto pit : points )
                    cgalpoints.emplace_back( pit[0], pit[1], pit[2] );
               return Surface( cgalpoints, num_neighbors );
          } ) , DOC(Surface,Surface,5 ))

          .def("set_proximity_ratio",     &Surface::set_proximity_ratio,     DOC(Surface,set_proximity_ratio ) ) 
          .def("set_smoothing_reduction", &Surface::set_smoothing_reduction, DOC(Surface,set_smoothing_reduction )  ) 
          .def("set_displacment_ratio",   &Surface::set_displacment_ratio,   DOC(Surface,set_displacment_ratio )  ) 
          
          .def("segment",[](Surface &self, double threshold)
          {
               auto face_markers = self.surface_segmentation(threshold);
               return py::array(py::cast(std::move(face_markers))); 
          },  py::return_value_policy::copy, py::arg("threshold")=60, DOC(Surface,segment))
  
          .def("get_faces", [](Surface &self){
               auto faces  = self.get_faces();
               return py::array(py::cast(std::move(faces))); 
          },   DOC(Surface,get_faces))

          .def("get_points", [](Surface &self){
               std::vector<Point_3>  cgalpoints = self.get_points();
               std::vector< std::array<double, 3>> points = wrapper_point_vector(cgalpoints);
               return py::array(py::cast(std::move(points)));
          },   DOC(Surface,get_points))
          
          .def("get_edges", [](Surface &self){
               auto edges  = self.get_edges();
               return py::array(py::cast(std::move(edges))); 
          } ,  DOC(Surface,get_edges))

          .def("__copy__", [](const Surface &self)
               { return Surface(self); })
          
          .def("copy",   [](const Surface &self)
               { return Surface(self); }, DOC(Surface,copy))
          
          .def("assign", &Surface::operator=)
        
          .def("keep_largest_connected_component", &Surface::keep_largest_connected_component, DOC(Surface, keep_largest_connected_component))
          .def("connected_components",  &Surface::connected_components, DOC(Surface, connected_components))
                
          .def("repair_self_intersections", &Surface::repair_self_intersections, 
               py::arg("volume_threshold") = 0.4, 
               py::arg("cap_threshold") = 170,
               py::arg("needle_threshold") = 1.5, 
               py::arg("collapse_threshold") = 0.3, DOC(Surface, repair_self_intersections))

          .def("remove_small_components", &Surface::remove_small_components, 
               py::arg("volume_threshold") = 30, DOC(Surface, remove_small_components))
        
          .def("implicit_surface", &Surface::implicit_surface<Domain,Surface_implicit_function>, 
               py::arg("implicit_function"), 
               py::arg("edge_length"),
               py::arg("bounding_radius"), 
               py::arg("error_bound") = 1.e-6, DOC(Surface, implicit_surface))

          .def("clip", py::overload_cast<double, double, double, double, bool>(&Surface::clip), 
               py::arg("x0"),
               py::arg("x1"),
               py::arg("x2"),
               py::arg("x3"),
               py::arg("preserve_manifold") = true, DOC(Surface, clip))
          
          .def("clip" , [](Surface &self, std::tuple<double,double,double,double> plane, bool preserve_manifold){
               self.clip( Wrapper_plane_3(plane), preserve_manifold);
          } )

          .def("clip" , [](Surface &self, std::tuple<double,double,double> point, std::tuple<double,double,double> vector, bool preserve_manifold){
               self.clip( Wrapper_point_3(point), Wrapper_vector_3(vector), preserve_manifold);
          } )

          .def("clip", py::overload_cast<Point_3, Vector_3, bool>(&Surface::clip), 
               py::arg("point"), 
               py::arg("vector"), 
               py::arg("preserve_manifold") = true, DOC(Surface, clip, 2))
        
          .def("clip", py::overload_cast<Plane_3, bool>(&Surface::clip), 
               py::arg("plane"), 
               py::arg("preserve_manifold") = true ,DOC(Surface, clip, 3) ) 
       
          .def("clip", py::overload_cast<Surface, bool, bool>(&Surface::clip), 
               py::arg("surface"), 
               py::arg("invert") = false, 
               py::arg("preserve_manifold") = true, DOC(Surface, clip, 4))
           
          .def("clip", py::overload_cast<Point_3, Vector_3, double, bool, bool>(&Surface::clip),
               py::arg("point"), 
               py::arg("vector"), 
               py::arg("radius"), 
               py::arg("invert") = false, 
               py::arg("preserve_manifold") = true, DOC(Surface, clip, 5))

          .def("get_slice", py::overload_cast<double, double, double, double>(&Surface::get_slice<Slice>), py::return_value_policy::move ,DOC(Surface, get_slice))
          .def("get_slice", py::overload_cast<Plane_3>(&Surface::get_slice<Slice>),py::return_value_policy::move,  DOC(Surface, get_slice, 2))
          .def("get_slice", [] ( Surface& self, std::tuple<double,double,double,double> plane){
               self.get_slice<Slice>( Wrapper_plane_3(plane));
          }, DOC(Surface,get_slice,3))
        
          .def("intersection",&Surface::surface_intersection,    DOC(Surface, surface_intersection))
          .def("union",       &Surface::surface_union,           DOC(Surface, surface_union))
          .def("difference",  &Surface::surface_difference,      DOC(Surface, surface_difference))

          .def("clear",  &Surface::clear,    DOC(Surface, clear))
          .def("span",   &Surface::span,     DOC(Surface, span)) 
          .def("save", []( Surface &self, std::string filename)
          {
               self.assert_non_empty_mesh();  
               self.set_outward_face_orientation();  
               std::string extension = filename.substr(filename.find_last_of(".")+1);
               std::ofstream out(filename);
               if( extension =="npz")
               { 
                    std::vector<std::array<double, 3>> points = wrapper_point_vector(self.get_points());
                    auto faces  = self.get_faces();
                    numpy_savez_wrapper(filename, points, faces);
               }
               else
                    self.save(filename);

          },py::arg("filename") , DOC(Surface, save) )    //&Surface::save,     DOC(Surface, save))


          //.def("save",   &Surface::save,     DOC(Surface, save))


          .def("fill_holes",            &Surface::fill_holes,         DOC(Surface, fill_holes))
          .def("triangulate_faces",     &Surface::triangulate_faces,  DOC(Surface, triangulate_faces))

          .def("isotropic_remeshing", py::overload_cast<double, unsigned int, bool>(&Surface::isotropic_remeshing),
               py::arg("edge_length")=1.0,
               py::arg("nb_iter")=5, 
               py::arg("protect_borders")=false, DOC(Surface, isotropic_remeshing))

          .def("isotropic_remeshing", py::overload_cast<double, unsigned int, double ,bool>(&Surface::isotropic_remeshing),
               py::arg("edge_length")=1.0,
               py::arg("nb_iter")=5, 
               py::arg("detect_border")=60,
               py::arg("protect_borders")=true, DOC(Surface, isotropic_remeshing,2))
          
          .def("isotropic_remeshing", py::overload_cast<double, unsigned int, std::vector<int>, bool>(&Surface::isotropic_remeshing),
               py::arg("edge_length")=1.0,
               py::arg("nb_iter")=5, 
               py::arg("face_markers"),
               py::arg("protect_borders")=true, DOC(Surface, isotropic_remeshing,2))
               
          .def("adjust_boundary",  py::overload_cast<double>(&Surface::adjust_boundary), DOC(Surface, adjust_boundary))
          .def("adjust_boundary",  py::overload_cast<Surface,double>( &Surface::adjust_boundary))        
          
          .def("smooth_laplacian", &Surface::smooth_laplacian,   DOC(Surface, smooth_laplacian))
          .def("smooth_taubin",    &Surface::smooth_taubin,      DOC(Surface, smooth_taubin))
          .def("smooth_shape",     &Surface::smooth_shape,       DOC(Surface, smooth_shape))
        
          .def("make_cube", py::overload_cast<double, double, double, double, double, double, double>(&Surface::make_cube),
               py::arg("x0"), 
               py::arg("y0"), 
               py::arg("z0"), 
               py::arg("x1"), 
               py::arg("y1"), 
               py::arg("z1"), 
               py::arg("edge_length"), DOC(Surface, make_cube, 2))

          .def("make_cube", py::overload_cast<Point_3, Point_3, double>(&Surface::make_cube), 
               py::arg("p1"), 
               py::arg("p2"),
               py::arg("edge_length"), DOC(Surface, make_cube))

          .def("make_cube", []( Surface& self, std::tuple<double,double,double> p1, std::tuple<double,double,double> p2, double edge_length){
               self.make_cube( Wrapper_point_3(p1), Wrapper_point_3(p2), edge_length);
          }, DOC(Surface,make_cube,3))
        
          .def("make_cone", py::overload_cast<double, double, double, double, double, double, double, double, double>(&Surface::make_cone), DOC(Surface, make_cone))
          .def("make_cone", py::overload_cast<Point_3, Point_3, double, double, double>(&Surface::make_cone), DOC(Surface, make_cone, 2))
          .def("make_cone", [](Surface &self, std::tuple<double,double,double> p0, std::tuple<double,double,double> p1, double r0, double r1,double edge_length) {
               self.make_cone( Wrapper_point_3(p0), Wrapper_point_3(p1), r0, r1, edge_length );
          })

          .def("make_cylinder", py::overload_cast<double, double, double, double, double, double, double, double>(&Surface::make_cylinder), DOC(Surface, make_cylinder))
          .def("make_cylinder", py::overload_cast<Point_3, Point_3, double, double>(&Surface::make_cylinder),                               DOC(Surface, make_cylinder, 2))
          .def("make_cylinder", [](Surface &self, std::tuple<double,double,double> p0, std::tuple<double,double,double> p1, double radius,double edge_length) {
               self.make_cone( Wrapper_point_3( p0), Wrapper_point_3(p1), radius, radius, edge_length );
          })

          .def("make_sphere", py::overload_cast< double, double, double, double, double, double >(&Surface::make_sphere<Domain>), 
               py::arg("x"), 
               py::arg("y"), 
               py::arg("z"), 
               py::arg("radius"), 
               py::arg("edge_length"), 
               py::arg("error_bound")=1.e-6, DOC(Surface, make_sphere))
        
          .def("make_sphere", py::overload_cast< Point_3, double, double, double>(&Surface::make_sphere<Domain>),
               py::arg("center"), 
               py::arg("radius"), 
               py::arg("edge_length"), 
               py::arg("error_bound")=1.e-6,   DOC(Surface, make_sphere, 2))

          .def( "make_sphere", [] ( Surface &self, std::tuple<double,double,double> center, double radius, double edge_length, double error_bound  ) {
               self.make_sphere<Domain>( Wrapper_point_3(center), radius, edge_length, error_bound );},  
               py::arg("center"), 
               py::arg("radius"), 
               py::arg("edge_length"), 
               py::arg("error_bound")=1.e-6)     

          .def("make_circle_in_plane", py::overload_cast<double, double, double, double, double, double, double, double>(&Surface::make_circle_in_plane),
               py::arg("x"), 
               py::arg("y"),
               py::arg("z"),
               py::arg("a"), 
               py::arg("b"),
               py::arg("c"),             
               py::arg("radius"), 
               py::arg("edge_length"), DOC(Surface,make_circle_in_plane,2) )
             
          .def("make_circle_in_plane", py::overload_cast<Point_3, Vector_3, double, double>(&Surface::make_circle_in_plane),
               py::arg("center"),
               py::arg("normal"),             
               py::arg("radius"), 
               py::arg("edge_length"), DOC(Surface,make_circle_in_plane) )             
             
          .def("make_circle_in_plane", [](Surface &self, std::tuple<double,double,double>  center, std::tuple<double,double,double> normal, double radius, double edge_length) {
               self.make_circle_in_plane(Wrapper_point_3( center ), Wrapper_vector_3(normal), radius , edge_length);
          }, py::arg("center"), py::arg("normal"), py::arg("radius"), py::arg("edge_length"), DOC(Surface,make_circle_in_plane,3) )             
             
          .def("is_point_inside", py::overload_cast<Point_3>(&Surface::is_point_inside), DOC(Surface, is_point_inside))
          .def("is_point_inside", [] (Surface &self,std::tuple<double,double,double> point ){
               self.is_point_inside( Wrapper_point_3(point)); 
          }, DOC(Surface, is_point_inside,2))

          .def("get_closest_points", py::overload_cast<Point_3,int>(&Surface::get_closest_points), 
               py::arg("point"), 
               py::arg("num") = 8, DOC(Surface, get_closest_points))

          .def("get_closest_points", [] ( Surface &self, std::tuple<double,double,double> point, int num){
               std::vector<Point_3> close_points = self.get_closest_points( Wrapper_point_3(point), num );
               return  py::array(py::cast(std::move(close_points)));}, py::arg("point"), py::arg("num") = 8,
               py::return_value_policy::copy, DOC(Surface, get_closest_points,2))

          .def("mean_curvature_flow", py::overload_cast<>(&Surface::mean_curvature_flow),            DOC(Surface, mean_curvature_flow))        
          .def("mean_curvature_flow", py::overload_cast<std::string>(&Surface::mean_curvature_flow), DOC(Surface, mean_curvature_flow,2)) //TODO
     
          .def("remove_degenerate_faces", &Surface::remove_degenerate_faces, DOC(Surface,remove_degenerate_faces) )  

          .def("get_shortest_surface_path", py::overload_cast<double, double, double, double, double, double>(&Surface::get_shortest_surface_path),
               DOC(Surface, get_shortest_surface_path))

          .def("get_shortest_surface_path", py::overload_cast<Point_3, Point_3>(&Surface::get_shortest_surface_path), DOC(Surface, get_shortest_surface_path, 2))
          .def("get_shortest_surface_path", [] ( Surface &self, std::tuple<double,double,double>  source , std::tuple<double,double,double>  target){
               std::vector<Point_3> cgalpath = self.get_shortest_surface_path( Wrapper_point_3(source) ,Wrapper_point_3(target) );
               std::vector<std::array<double,3>> path = wrapper_point_vector(cgalpath);
               return  py::array(py::cast(std::move(path)));}, py::return_value_policy::copy, DOC(Surface, get_shortest_surface_path, 3))

          .def("embed", &Surface::embed, py::arg("other"), 
               py::arg("adjustment") = .4, 
               py::arg("smoothing") = 0.4, 
               py::arg("max_iter") = 50, DOC(Surface, embed))
        
          .def("enclose", &Surface::enclose, 
               py::arg("other"), 
               py::arg("adjustment") = .4, 
               py::arg("smoothing") = 0.4, 
               py::arg("max_iter") = 50, DOC(Surface, enclose))
        
          .def("expose", &Surface::expose, 
               py::arg("other"), 
               py::arg("adjustment") = .4, 
               py::arg("smoothing") = 0.4, 
               py::arg("max_iter") = 50, DOC(Surface, expose))
        
          .def("separate", &Surface::separate,
               py::arg("other"), 
               py::arg("adjustment"), 
               py::arg("smoothing") = 0.4, 
               py::arg("max_iter") = 50, DOC(Surface, separate))

          .def("collapse_edges", py::overload_cast<const double>(&Surface::collapse_edges), DOC(Surface, collapse_edges))
          .def("collapse_edges", py::overload_cast<>(&Surface::collapse_edges), DOC(Surface, collapse_edges, 2))
          
          .def("split_edges", &Surface::split_edges, DOC(Surface, split_edges))
        
          .def("get_perpendicular_cut", py::overload_cast<Plane_3, double>(&Surface::get_perpendicular_cut),
               py::arg("plane") , 
               py::arg("radius")=0.0, py::return_value_policy::move, DOC(Surface,get_perpendicular_cut,2)) 
        
          .def("get_perpendicular_cut", py::overload_cast<Point_3, double>(&Surface::get_perpendicular_cut), 
               py::arg("query"), 
               py::arg("radius")=0.0, py::return_value_policy::move, DOC(Surface,get_perpendicular_cut)) 

          .def("get_perpendicular_cut",[]( Surface &self,std::tuple<double,double,double> query, double radius ){ 
               self.get_perpendicular_cut( Wrapper_point_3(query), radius );
          })

          .def("extension", py::overload_cast<double, double, double, double, double, double, bool>(&Surface::cylindrical_extension), 
               py::arg("x"), 
               py::arg("y"), 
               py::arg("z"),
               py::arg("radius"),
               py::arg("length"),
               py::arg("edge_length"),
               py::arg("use_normal"), pybind11::return_value_policy::move, DOC(Surface, cylindrical_extension)) //TODO test return policy
             
          .def("extension", py::overload_cast<const Point_3 &, double, double, double, bool>(&Surface::cylindrical_extension), py::arg("point"),
               py::arg("radius"),
               py::arg("length"),
               py::arg("edge_length"),
               py::arg("use_normal"),  pybind11::return_value_policy::move, DOC(Surface, cylindrical_extension, 2))
          
          .def("extension", [](Surface &self, std::tuple<double,double,double> point, double radius, double length, double edge_length, bool use_normal ){
               return self.cylindrical_extension(Wrapper_point_3(point), radius, length, edge_length,use_normal);
          })     
        
          .def("connection", &Surface::cylindrical_connection, DOC(Surface, cylindrical_connection))

          .def("separate_narrow_gaps", &Surface::separate_narrow_gaps, 
               py::arg("adjustment") = 0.5, 
               py::arg("smoothing") = 0.5, 
               py::arg("max_iter") = 100, DOC(Surface, separate_narrow_gaps))
             
          .def("separate_close_vertices", &Surface::separate_close_vertices, 
               py::arg("adjustment") = 0.6, 
               py::arg("max_iter") = 100, DOC(Surface, separate_close_vertices))

          .def("convex_hull", &Surface::convex_hull, DOC(Surface, convex_hull))

          .def("num_faces",        &Surface::num_faces,     DOC(Surface, num_faces))
          .def("num_edges",        &Surface::num_edges,     DOC(Surface, num_edges))
          .def("num_vertices",     &Surface::num_vertices,  DOC(Surface, num_vertices))
          .def("num_self_intersections", &Surface::num_self_intersections, DOC(Surface, num_self_intersections))

          .def("distance", py::overload_cast<Point_3>(&Surface::distance_to_point), DOC(Surface, distance_to_point))
          .def("distance",[] ( Surface &self, std::tuple<double,double,double> point ){
               return self.distance_to_point(Wrapper_point_3(point));
          },  DOC(Surface, distance_to_point,2) )
          
          .def("centeroid", &Surface::centeroid,         DOC(Surface, centeroid)) // TODO 
          
          .def("area",      &Surface::area,              DOC(Surface, area))
          .def("volume",    &Surface::volume,            DOC(Surface, volume));

     py::class_<Domain, std::shared_ptr<Domain>>(m, "Domain", DOC(Domain))
          .def(py::init<Surface &, double>(), 
               py::arg("surface"), 
               py::arg("error_bound") = 1.e-3, DOC(Domain, Domain))
             
          .def(py::init<std::vector<Surface>, double>(), 
               py::arg("surfaces"), 
               py::arg("error_bound") = 1.e-3, DOC(Domain, Domain, 2))
        
          .def(py::init<std::vector<Surface>, std::shared_ptr<AbstractMap>, double>(), 
             py::arg("surfaces"), 
             py::arg("map"), 
             py::arg("error_bound") = 1.e-3, DOC(Domain, Domain, 3))
        
          .def(py::init<std::string &,double , Surface>(), 
             py::arg("filename"), 
             py::arg("error_bound") = 1.e-3, 
             py::arg("surface")= Surface(), DOC(Domain, Domain, 4))         
        
          .def("tetrahedral_remeshing",&Domain::tetrahedral_remeshing,
             py::arg("edge_length")=1.0,
             py::arg("nb_iter")=5, 
             py::arg("protect_borders")=false, DOC(Domain,tetrahedral_remeshing))       

          .def("get_facets", []( Domain &self, bool exclude_unmarked, int first_index){  
             auto facets  = self.get_facets(exclude_unmarked, first_index);
             return py::array(py::cast(std::move(facets)));
          }, py::arg("exclude_unmarked")=true , py::arg("first_index")=1, py::return_value_policy::copy, DOC(Domain,get_facets) ) 
       
          .def("get_facet_tags", []( Domain &self,  bool exclude_unmarked){ 
               auto facet_tags = self.get_facet_tags(exclude_unmarked);
               return py::array(py::cast(std::move(facet_tags)));
          }, py::arg("exclude_unmarked")=true  , py::return_value_policy::copy, DOC(Domain,get_facet_tags) )      
        
          .def("get_cell_tags", []( Domain &self){ 
               auto cell_tags = self.get_cell_tags();
               return py::array(py::cast(std::move(cell_tags)));
          }, py::return_value_policy::copy, DOC(Domain,get_cell_tags) )   
        
          .def("get_cells", []( Domain &self, int first_index){
               auto cells = self.get_cells(first_index);
               return  py::array(py::cast(std::move(cells)));       
          }, py::arg("first_index")=1, py::return_value_policy::copy,DOC(Domain,get_cells) ) 

          .def("get_points",[]( Domain &self){
               auto points = self.get_points();
               return  py::array(py::cast(std::move(points)));       
          }, py::return_value_policy::copy, DOC(Domain,get_points) ) 

          .def("get_borders",  &Domain::get_borders , DOC(Domain,get_borders)  )  
          .def("get_features", &Domain::get_features, DOC(Domain,get_features) )        

          .def("create_mesh", py::overload_cast<double, double, double, double, double, double, double, double>(&Domain::create_mesh),
               py::arg("edge_size"), 
               py::arg("cell_size"), 
               py::arg("facet_size"),
               py::arg("facet_angle"), 
               py::arg("facet_distance"), 
               py::arg("cell_radius_edge_ratio"),
               py::arg("min_edge_size"),
               py::arg("min_cell_size"), DOC(Domain, create_mesh))

          .def("create_mesh", py::overload_cast<double>(&Domain::create_mesh),  DOC(Domain, create_mesh, 2))
          .def("create_mesh", py::overload_cast<>(&Domain::create_mesh),        DOC(Domain, create_mesh, 3))
        
          .def("radius_ratios_min_max",   &Domain::radius_ratios_min_max,   DOC(Domain, radius_ratios_min_max))
          .def("dihedral_angles_min_max", &Domain::dihedral_angles_min_max, DOC(Domain, dihedral_angles_min_max))

          .def("radius_ratios",   &Domain::radius_ratios,   DOC(Domain, radius_ratios))
          .def("dihedral_angles", &Domain::dihedral_angles, DOC(Domain, dihedral_angles))

          .def("get_boundary", &Domain::get_boundary, 
               py::arg("tag") = 0, DOC(Domain, get_boundary))
             
          .def("get_boundaries", &Domain::get_boundaries, DOC(Domain, get_boundaries))

          .def("get_borders",    &Domain::get_borders,            DOC(Domain, get_borders)) 
          .def("get_interface",  &Domain::get_interface, DOC(Domain, get_interface))
          .def("get_curve_tags", &Domain::get_curve_tags,         DOC(Domain, get_curve_tags))
          .def("get_patches",    &Domain::get_patches,            DOC(Domain, get_patches))
          .def("get_subdomains", &Domain::get_subdomains,         DOC(Domain, get_subdomains))
        
          .def("validate_mesh", &Domain::validate_mesh, DOC(Domain, check_mesh_connections))
        
          .def("init_triangulation", py::overload_cast<Surface>(&Domain::init_triangulation),              DOC(Domain,init_triangulation)  )           
          .def("init_triangulation", py::overload_cast<std::vector<Surface>>(&Domain::init_triangulation), DOC(Domain,init_triangulation,2) ) 
       
          .def("refine", &Domain::refine, DOC(Domain,refine))   
     
          .def("exude", &Domain::exude, 
               py::arg("time_limit") = 10., 
               py::arg("sliver_bound") = 0, DOC(Domain, exude))
        
          .def("perturb", &Domain::perturb, 
               py::arg("time_limit") = 10., 
                py::arg("sliver_bound") = 0, DOC(Domain, perturb))
         
          .def("lloyd", &Domain::lloyd, 
               py::arg("time_limit") = 10.,
               py::arg("max_iter") = 0,
               py::arg("convergence") = 0.02,
               py::arg("freeze_bound") = 0.01,
               py::arg("do_freeze") = true, DOC(Domain, lloyd))

         .def("odt", &Domain::odt, 
               py::arg("time_limit") = 10.,
               py::arg("max_iter") = 0,
               py::arg("convergence") = 0.01,
               py::arg("freeze_bound") = 0.001,
               py::arg("do_freeze") = false, DOC(Domain, odt))

          .def("add_sharp_border_edges", py::overload_cast<Surface &, double>(&Domain::add_sharp_border_edges), 
               py::arg("surface"),
               py::arg("threshold") = 60, DOC(Domain, add_sharp_border_edges))
         
          .def("add_sharp_border_edges", py::overload_cast<Surface &, Plane_3, double>(&Domain::add_sharp_border_edges), 
               py::arg("surface"),
               py::arg("plane"),
               py::arg("threshold") = 60, DOC(Domain, add_sharp_border_edges, 2))

          .def("add_sharp_border_edges", []( Domain &self, Surface &surface, std::tuple<double,double,double,double> plane, double threshold){
               self.add_sharp_border_edges(surface, Wrapper_plane_3(plane), threshold);
          } , py::arg("surface"), py::arg("plane"), py::arg("threshold") = 60, DOC(Domain, add_sharp_border_edges, 4))

        
         .def("add_sharp_border_edges", py::overload_cast<Surface &, Surface&, double>(&Domain::add_sharp_border_edges), 
               py::arg("surface"),
               py::arg("clip"),
               py::arg("threshold") = 60, DOC(Domain, add_sharp_border_edges, 3))           
     
          .def("add_sharp_border_edges", []( Domain &self, Surface &surface, std::tuple<double,double,double, double> plane, double threshold  ){
               self.add_sharp_border_edges( surface, Wrapper_plane_3(plane), threshold);
          })      

         .def("clear_borders",     &Domain::clear_borders,  DOC(Domain, clear_borders))
         .def("clear_features",    &Domain::clear_features, DOC(Domain, clear_features))
        
         .def("remove_subdomain", py::overload_cast<std::vector<int>>(&Domain::remove_subdomain),   DOC(Domain, remove_subdomain))
         .def("remove_subdomain", py::overload_cast<int>(&Domain::remove_subdomain),                DOC(Domain, remove_subdomain, 2))

         .def("num_cells",      &Domain::num_cells,      DOC(Domain, num_cells))
         .def("num_subdomains", &Domain::num_subdomains, DOC(Domain, num_subdomains))
         .def("num_curves",     &Domain::num_curves,     DOC(Domain, num_curves))
         .def("num_patches",    &Domain::num_patches,    DOC(Domain, num_patches))
         .def("num_surfaces",   &Domain::num_surfaces,   DOC(Domain, num_surfaces))
         .def("num_facets",     &Domain::num_facets,     DOC(Domain, num_facets))
         .def("num_vertices",   &Domain::num_vertices,   DOC(Domain, num_vertices))
        
         //.def("write_facet_data", &Domain::write_facet_data, DOC(Domain, write_facet_data)) 
        
         .def("get_collision_distances", &Domain::get_collision_distances,  DOC(Domain, get_collision_distances))
         .def("get_collision_spheres"  , &Domain::get_collision_spheres,    DOC(Domain, get_collision_spheres))
              
          .def("boundary_segmentations", py::overload_cast<std::pair<int, int>, double>(&Domain::boundary_segmentations),
               py::arg("interface"),
               py::arg("angle_in_degree") = 85, DOC(Domain, boundary_segmentations))

          .def("boundary_segmentations", py::overload_cast<int, double>(&Domain::boundary_segmentations),
               py::arg("subdomain_tag"),
               py::arg("angle_in_degree") = 85, DOC(Domain, boundary_segmentations, 2))

          .def("boundary_segmentations", py::overload_cast<double>(&Domain::boundary_segmentations),
               py::arg("angle_in_degree") = 85, DOC(Domain, boundary_segmentations, 3))

          .def("add_feature", &Domain::add_feature, DOC(Domain, add_feature))
          .def("add_border", &Domain::add_border,   DOC(Domain, add_border))
        
          // .def("save", py::overload_cast<std::string, bool ,bool>(&Domain::save),
          //      py::arg("OutPath"),
          //      py::arg("exclude_unmarked_facets") = false,
          //      py::arg("save_edge_features") = false, DOC(Domain, save));

          .def("save", []( Domain &self, std::string filename, bool exclude_unmarked, bool save_edge_features)
          {
               self.assert_non_empty_mesh_object();
               std::string extension = filename.substr(filename.find_last_of(".")+1);
               std::ofstream out(filename);
               if( extension =="npz")
               { 
                    auto points = self.get_points(); 
                    auto cells  = self.get_cells(0); // first index
                    auto cell_tags = self.get_cell_tags();
                    auto facets = self.get_facets(exclude_unmarked,0);
                    auto facet_tags = self.get_facet_tags(exclude_unmarked);
                    if( save_edge_features)
                    {
                         auto edges = self.get_edges(exclude_unmarked,0 );
                         auto edge_tags = self.get_edge_tags(exclude_unmarked);
                         numpy_savez_wrapper(filename, points, cells, cell_tags,
                              facets, facet_tags, edges, edge_tags);
                    }
                    else 
                         numpy_savez_wrapper(filename, points, cells, cell_tags, facets, facet_tags);               
               }
               else
                    self.save(filename, exclude_unmarked, save_edge_features);



          } , py::arg("filename"), py::arg("exclude_unmarked") = false, py::arg("save_edge_features") = false, DOC(Domain, save) );
          // TODO GENERAL CHECK first_index


     m.def("convex_hull", &Wrapper_convex_hull, DOC(convex_hull)); // TODO

     m.def("smooth_polyline", &smooth_polyline<Point_2,Vector_2>,     DOC(smooth_polylines));
     m.def("smooth_polyline", &smooth_polyline<Point_3,Vector_3>,     DOC(smooth_polylines,2)); 
     m.def("inside_polygon",  &inside_polygon<Point_2>,               DOC(inside_polygon));  


     m.def("enclose", py::overload_cast<Surface &, Surface &, double, double, int>(&enclose),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("adjustment") = .4,
          py::arg("smoothing") = 0.4,
          py::arg("max_iter") = 50, DOC(enclose)); 

     m.def("embed", py::overload_cast<Surface &, Surface &, double, double, int>(&embed),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("adjustment") = -.4,
          py::arg("smoothing") = 0.4,
          py::arg("max_iter") = 50, DOC(embed)); 

     m.def("expose", py::overload_cast<Surface &, Surface &, double, double, int>(&expose),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("adjustment") = -.4,
          py::arg("smoothing") = 0.4,
          py::arg("max_iter") = 50, DOC(expose)); 
    
     m.def("separate_overlapping_surfaces", py::overload_cast<Surface &, Surface &, Surface &, double, double, int>(&separate_surface_overlapp),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("other"),
          py::arg("edge_movement") = -.4,
          py::arg("smoothing") = 0.4,
          py::arg("max_iter") = 50, DOC(separate_surface_overlapp));

     m.def("separate_overlapping_surfaces", py::overload_cast<Surface &, Surface &, double, double, int>(&separate_surface_overlapp),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("edge_movement") = -.4,
          py::arg("smoothing") = 0.4,
          py::arg("max_iter") = 50, DOC(separate_surface_overlapp, 2));


     m.def("separate_close_surfaces", py::overload_cast<Surface &, Surface &, Surface &, double, double, int>(&separate_close_surfaces),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("other"),
          py::arg("edge_movement") = -.4,
          py::arg("smoothing") = 0.4,
          py::arg("max_iter") = 50, DOC(separate_close_surfaces));

     m.def("separate_close_surfaces", py::overload_cast<Surface&, Surface&, double, double, int>(&separate_close_surfaces),
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("edge_movement") = -.4,
          py::arg("smoothing") = 0.4,
          py::arg("max_iter") = 50, DOC(separate_close_surfaces, 2));

     m.def("union_partially_overlapping_surfaces", &union_partially_overlapping_surfaces,
          py::arg("surf1"),
          py::arg("surf2"),
          py::arg("angle_in_degress") = 36.87,
          py::arg("adjustment") = 2.0,
          py::arg("smoothing") = 0.2,
          py::arg("max_iter") = 50, DOC(union_partially_overlapping_surfaces)); 
}
