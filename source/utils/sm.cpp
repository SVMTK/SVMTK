#include <vector>
#include <fstream>
#include <limits>
#include <boost/foreach.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <pybind11/pybind11.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;

namespace py = pybind11;


struct Surface_mesh {
    Polyhedron poly;

    Surface_mesh(const std::string &filename) {
        std::ifstream input(filename);
        if (!input || !(input >> poly) || poly.empty()
                   || !CGAL::is_triangle_mesh(poly))
        {
            std::cerr << "Not a valid input file." << std::endl;
            // TODO: Raise error
        }
    }

    void save(const std::string &filename) {
        /* std::ofstream out(output_filename.c_str()); */
        std::ofstream out(filename);
        out << poly;
    }
};


PYBIND11_MODULE(sm, m) {
    py::class_<Surface_mesh>(m, "Surface_mesh")
        .def(py::init<const std::string &>())
        .def("save", &Surface_mesh::save);
        
}
