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

//#include "read_polygons_STL.h"


namespace utils {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

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
        return 1;
    }


    template<typename Mesh>     // seperate because -> typedef mesh in CGALSurface
    bool load_surface(Mesh& mesh, std::string filename) {
        std::ifstream input(filename);

        // This seems unnecessary
        std::string file(filename);
        std::string extension = file.substr(file.find_last_of(".") + 1);

        if (!input) {
            std::cerr << "Cannot open file " << std::endl;
            return false;
        }

        typedef typename Mesh::Point Point_3;
        std::vector<Point_3> points;
        std::vector< std::vector<std::size_t> > polygons;

        if (extension == "off") {
            std::cout << "reading off" << std::endl;
            if (!CGAL::read_OFF(input, points, polygons)) {
                std::cerr << "Error parsing the OFF file " << std::endl;
                return false;
            }
        }
        else if (extension == "stl") {
            // if (!read_polygons_STL(input, points, polygons)) { // TODO: the connection causes errors 
            //     std::cerr << "Error parsing the STL file " << std::endl;
            //     return false;
            // }
        }
        else {
            std::cerr << "Error unkown file extension" << std::endl;
            return false;
        }

        CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,mesh);

        if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh))) {
            CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
        }
        return true;
    }

}

#endif
