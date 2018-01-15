#include "surface_mesh.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;


PYBIND11_MODULE(brainmesh, m)  {
    py::class_<Surface_mesh>(m, "Surface_mesh")
        .def(py::init<const std::string &, bool>())
        .def("save", &Surface_mesh::save)
        .def("stitch", &Surface_mesh::stitch)
        .def("num_selfintersections", &Surface_mesh::num_selfintersections)
        .def("remesh", &Surface_mesh::remesh)
        .def("fill_holes", &Surface_mesh::fill_holes);
}
