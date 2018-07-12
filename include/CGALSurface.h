#ifndef CGALSurface_H
#define CGALSurface_H

#define BOOST_PARAMETER_MAX_ARITY 12

#include "utils.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

// Isotropic remeshing
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

// Edge collapse
// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h> 

// copy_face_graph
#include <CGAL/boost/graph/copy_face_graph.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel; // Change to exact -> copy face_ graph to inexact ??
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef Kernel::Vector_3 Vector_3;
typedef boost::graph_traits<Mesh>::face_descriptor        face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor      vertex_descriptor;


class CGALSurface {
    /* typedef std::vector<vertex_descriptor>                    vertex_vector; */
    /* typedef CGAL::Side_of_triangle_mesh<Mesh,Kernel> Inside; // */
    /*     typedef CGAL::Surface_mesh_default_triangulation_3 Tr; */
    /*     typedef Tr::Geom_traits GT; */

    private:
        Mesh mesh;

    public:
        CGALSurface(const std::string f);

        template<typename Polyhedron_3>
        void get_polyhedron(Polyhedron_3 &polyhedron_3);

    /* CGALSurface( Implicit_function implicit_function, */
    /*              double bounding_sphere_radius, */
    /*              double angular_bound, */
    /*              double radius_bound, */
    /*              double distance_bound); */


    /*     void operator^=( CGALSurface& other); */

    /*     void operator+=( CGALSurface& other ); */

    /*     void operator-=( CGALSurface& other ); */

        Mesh& get_mesh();

        unsigned int fill_holes();

        bool triangulate_faces();

        void stitch_borders(); 

    /*     void insert_surface(CGALSurface& surface); */    // TODO

        void isotropic_remeshing(const double, const unsigned int, const bool);

        void adjust_boundary(const double);

        void smooth_laplacian(const double);

        template<typename InputIterator>
        void adjusting_boundary_region(InputIterator, InputIterator, const double);

        template<typename InputIterator>
        void smooth_laplacian_region(InputIterator, InputIterator, const int);

    /*     vertex_vector points_inside(CGALSurface& other); */

    /*     vertex_vector points_outside(CGALSurface& other); */

    /*     Mesh& get_mesh(); */
    /*     //Polyhedron& polyhedron(); */

        bool self_intersections();

        void save(const std::string);

        int collapse_edges(const double stop_ratio);

        void preprocess(const double, const int);

        /* void fair(CGALSurface::vertex_vector vector); */

    /*     template<typename Implicit_function> */
    /*     CGALSurface(Implicit_function implicit_function, */
    /*          double bounding_sphere_radius, */
    /*          double angular_bound=30., */
    /*          double radius_bound=0.1 , */
    /*          double distance_bound=0.1   ); */

    /*     ~CGALSurface(){} */
};

CGALSurface::CGALSurface(const std::string f) {
    utils::read_off(mesh, f);
}


template<typename Polyhedron_3>
void CGALSurface::get_polyhedron(Polyhedron_3 &polyhedron_3) {
    CGAL::copy_face_graph(mesh, polyhedron_3);
}


Mesh& CGALSurface::get_mesh() {
    return mesh;
}

unsigned int CGALSurface::fill_holes() {
    unsigned int nb_holes = 0;
    for (auto h: halfedges(mesh)) {
        if(is_border(h, mesh)) {
            std::vector<face_descriptor> patch_facets;
            std::vector<vertex_descriptor> patch_vertices;
            bool success = CGAL::cpp11::get<0>(
            CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(mesh, h, std::back_inserter(patch_facets),
                std::back_inserter(patch_vertices),
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).geom_traits(Kernel())) );
            std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
            std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
            std::cout << "  Is fairing successful: " << success << std::endl;
            nb_holes++;
        }
    }
    std::cout << std::endl;
    std::cout << nb_holes << " holes have been filled" << std::endl;
    return nb_holes;
}


bool CGALSurface::triangulate_faces() {
    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    for (auto fit: faces(mesh)) {
        if (next(next(halfedge(fit,  mesh), mesh), mesh) != prev(halfedge(fit, mesh), mesh)) {
            std::cerr << "Error: non-triangular face left in mesh." << std::endl;
        }
    }
    return true;
}


void CGALSurface::stitch_borders() {
    CGAL::Polygon_mesh_processing::stitch_borders(mesh);
}


void CGALSurface::isotropic_remeshing(
    const double target_edge_length,
    const unsigned int nb_iter,
    const bool protect_border) {
    CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length, mesh);
    CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh),
            target_edge_length, mesh,
            CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
            .protect_constraints(protect_border));
}


void CGALSurface::adjust_boundary(const double c) {
    Mesh::Vertex_range::iterator vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    CGALSurface::adjusting_boundary_region(vb, ve, c);
}


void CGALSurface::smooth_laplacian(const double c) {
    Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    CGALSurface::smooth_laplacian_region(vb, ve, c);
}


template<typename InputIterator>
void CGALSurface::adjusting_boundary_region(
        InputIterator begin,
        InputIterator end,
        const double c) {
    std::vector<std::pair<vertex_descriptor, Point_3> > smoothed; //rename
    for ( ; begin != end; ++begin) {
        Vector_3 delta= CGAL::Polygon_mesh_processing::compute_vertex_normal(*begin,mesh);
        Point_3 p = mesh.point(*begin) + c*delta;
        smoothed.push_back(std::make_pair(*begin, p));
    }
    for (const std::pair<vertex_descriptor, Point_3> s: smoothed) {
        mesh.point(s.first) = s.second;
    }
}


template<typename InputIterator>
void CGALSurface::smooth_laplacian_region(InputIterator begin, InputIterator end, const int c) {
    std::vector<std::pair<vertex_descriptor, Point_3> > smoothed;
    for ( ; begin != end; ++begin) {
        Point_3 current = mesh.point(*begin);
        Vector_3 delta = CGAL::NULL_VECTOR;
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(*begin),mesh), done(vbegin);
        do {
            delta += Vector_3(mesh.point(*vbegin) - current);
            *vbegin++;
        } while(vbegin != done);
        Point_3 p = current + c*delta/mesh.degree(*begin);
        smoothed.push_back(std::make_pair(*begin, p));
    }
    for (const std::pair<vertex_descriptor, Point_3> s : smoothed) {
        mesh.point(s.first) = s.second;
    }
}

bool CGALSurface::self_intersections() {
    return CGAL::Polygon_mesh_processing::does_self_intersect(mesh,
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(
                    get(CGAL::vertex_point, mesh)));
}


void CGALSurface::save(const std::string outpath) {
    utils::save_off(mesh, outpath);
};


int CGALSurface::collapse_edges(const double stop_ratio) {
    namespace SMS = CGAL::Surface_mesh_simplification;
    SMS::Count_ratio_stop_predicate<Mesh> stop(stop_ratio);

    const int r = SMS::edge_collapse(
    mesh,
    stop,
    CGAL::parameters::get_cost(SMS::Edge_length_cost <Mesh>())
        .get_placement(SMS::Midpoint_placement<Mesh>()));
        /* .visitor(vis)); */
    return r;
}

void CGALSurface::preprocess(const double target_edge_length, const int nb_iter) {
    CGALSurface::triangulate_faces();
    CGALSurface::isotropic_remeshing(target_edge_length, nb_iter, false);
};

#endif
