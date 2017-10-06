#include <fstream>
#include <vector>

#include "fixpolyhedron.h"
#include "remeshing.h"
#include "edge_collapse.h"
#include "hole_filling.h"
#include "simplify.h"

#include "stitching.h"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc, char* argv[])
{
  const char* filename = argv[1];
  const char* outputfile = argv[2];

  Mesh mesh;
  if (fix_polyhedron<Mesh, K>(filename, mesh)) {
    double target_edge_length = 0.04;
    unsigned int nb_iter = 3;
    edge_collapse(mesh);
    fill_holes(mesh);
    simplify(mesh);
    stitching(mesh);
    remeshing(mesh, target_edge_length, nb_iter);

    std::ofstream out(outputfile);
    out << mesh;
    out.close();
    return 0;
  }
  return 1;
}
