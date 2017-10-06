#include <fstream>
#include <vector>

#include "fixpolyhedron.h"
#include "remeshing.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/pig.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (fix_polyhedron<Mesh, K>(filename, mesh)) {
    // if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    //   std::cerr << "Not a valid input file." << std::endl;
    //   return 1;
    // }

    double target_edge_length = 0.04;
    unsigned int nb_iter = 3;

    remeshing(mesh, target_edge_length, nb_iter);
    return 0;
  }
  return 1;
}
