#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = argv[1];
  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }
  bool intersecting = PMP::does_self_intersect(mesh,
      PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));
  std::cout
    << (intersecting ? "There are self-intersections." : "There is no self-intersection.")
    << std::endl;
  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  PMP::self_intersections(mesh, std::back_inserter(intersected_tris));
  std::cout << intersected_tris.size() << " pairs of triangles intersect." << std::endl;

  std::vector<std::pair<face_descriptor, face_descriptor>>::iterator it;
  for(it = intersected_tris.begin(); it != intersected_tris.end(); it++) {
      std::cout << it->first << std::endl;
      std::cout << it->second << std::endl;
      std::cout << std::endl;
    // std::cout << it[0] << std::endl;
  }
  
  return 0;
}
