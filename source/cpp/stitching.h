#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

int stitching(Mesh& mesh) {
    std::cout << "Before stitching : " << std::endl;
    std::cout << "\t Number of vertices  :\t" << mesh.number_of_vertices() << std::endl;
    std::cout << "\t Number of halfedges :\t" << mesh.number_of_halfedges() << std::endl;
    std::cout << "\t Number of facets    :\t" << mesh.number_of_faces() << std::endl;

    CGAL::Polygon_mesh_processing::stitch_borders(mesh);
    std::cout << "Stitching done : " << std::endl;
    std::cout << "\t Number of vertices  :\t" << mesh.number_of_vertices() << std::endl;
    std::cout << "\t Number of halfedges :\t" << mesh.number_of_halfedges() << std::endl;
    std::cout << "\t Number of facets    :\t" << mesh.number_of_faces() << std::endl;

    return 0;
}

