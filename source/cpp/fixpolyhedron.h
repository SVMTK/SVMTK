#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <vector>
#include <fstream>
#include <iostream>

// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL::Polyhedron_3<K> Polyhedron;

template< typename Polyhedron , typename Kernel >
bool fix_polyhedron(const char* filename, Polyhedron& polyhedron) {
    typedef typename Kernel::Point_3 Point;
    std::ifstream input(filename);
  
    if (!input) {
        std::cerr << "Cannot open file " << std::endl;
        return false;
    }
    
    std::vector<Point> points;
    std::vector< std::vector<std::size_t> > polygons;
    if (!CGAL::read_OFF(input, points, polygons)) {
        std::cerr << "Error parsing the OFF file " << std::endl;
        return false;
    }
    
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
        points,
        polygons,
        polyhedron
    );
  
    if (CGAL::is_closed(polyhedron) && 
        (!CGAL::Polygon_mesh_processing::is_outward_oriented(polyhedron))) {
        std::cout<< "reverse_face_orientation"<< std::endl;
        CGAL::Polygon_mesh_processing::reverse_face_orientations(polyhedron); 
    }
    
    return true;
}














/*
int fixpolyhedron()
{
  const char* filename = (argc > 1) ? argv[1] : "data/tet-shuffled.off";
  std::ifstream input(filename);
  if (!input)
  {
    std::cerr << "Cannot open file " << std::endl;
    return 1;
  }
  std::vector<K::Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;
  if (!CGAL::read_OFF(input, points, polygons))
  {
    std::cerr << "Error parsing the OFF file " << std::endl;
    return 1;
  }
  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
  Polyhedron mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);
  if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
    CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
  std::ofstream out("tet-oriented1.off");
  out << mesh;
  out.close();
  CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
  std::ofstream out2("tet-oriented2.off");
  out2 << mesh;
  out2.close();
  return 0;
}
*/
