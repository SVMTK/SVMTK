#ifndef __CGAL_UTILS_H
#define __CGAL_UTILS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <vector>
#include <fstream>
#include <iostream>

namespace utils {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
    /* typedef CGAL::Polyhedron_3<Kernel> Polyhedron; */


    template<typename Mesh>
    int read_off(Mesh &mesh, const std::string filename) {
        std::ifstream input(filename);
        if (!input) {
            std::cout << "Cannot open file " << std::endl;
            return 1;
        }

        std::vector<Kernel::Point_3> points;
        std::vector< std::vector<std::size_t> > polygons;
        if (!CGAL::read_OFF(input, points, polygons)) {
            std::cout << "Error parsing the OFF file " << std::endl;
            return 1;
        }
        CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);
        if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh))) {
            CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
        }
        return 0;
    }


    template<typename Mesh>
    int save_off(const Mesh &mesh, const std::string outpath) {
        std::ofstream out(outpath);
        out << mesh;
        out.close();
    }

}

#endif
