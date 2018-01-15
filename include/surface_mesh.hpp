#include <fstream>
#include <iostream>
#include <vector>
#include <limits>

#include "edge_collapse.hpp"

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>

#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor       face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor      vertex_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;


struct halfedge2edge {
    halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges) :
            m_mesh(m), m_edges(edges) {}

    void operator()(const halfedge_descriptor& h) const {
        m_edges.push_back(edge(h, m_mesh));
    }

    const Mesh& m_mesh;
    std::vector<edge_descriptor>& m_edges;
};


int triangulate_faces(Mesh& mesh) {
    PMP::triangulate_faces(mesh);
  
    // Confirm that all faces are triangles.
    BOOST_FOREACH(boost::graph_traits<Mesh>::face_descriptor fit, faces(mesh)) {
        if (next(next(halfedge(fit, mesh), mesh), mesh) !=  prev(halfedge(fit, mesh), mesh)) {
            std::cerr << "Error: non-triangular face left in mesh." << std::endl;
        }
    }
  
    return 0;
}

struct Surface_mesh {
    Mesh surface_mesh;

    Surface_mesh(const std::string &filename, bool rearrange) {
        std::ifstream input(filename);
        if (rearrange) {
            if (!input || !(input >> surface_mesh))
            {
                std::cerr << "Not a valid input file." << std::endl;
                // TODO: Raise error
            }
        } else {
            std::vector<K::Point_3> points;
            std::vector< std::vector<std::size_t> > polygons;
            if (!CGAL::read_OFF(input, points, polygons)) {
                std::cerr << "Error parsing the OFF file " << std::endl;
            }
            
            PMP::orient_polygon_soup(points, polygons);
            PMP::polygon_soup_to_polygon_mesh(
                points,
                polygons,
                surface_mesh 
            );
      
            if (CGAL::is_closed(surface_mesh) && 
                (!PMP::is_outward_oriented(surface_mesh))) {
                std::cout<< "reverse_face_orientation"<< std::endl;
                PMP::reverse_face_orientations(surface_mesh); 
            }
            edge_collapse(surface_mesh);
        }
    }

    void save(const std::string &filename) {
        /* std::ofstream out(output_filename.c_str()); */
        std::ofstream out(filename);
        out << surface_mesh;
    }

    void stitch() {
        std::cout << "Before stitching : " << std::endl;
        std::cout << "\t Number of vertices  :\t" << surface_mesh.number_of_vertices() << std::endl;
        std::cout << "\t Number of halfedges :\t" << surface_mesh.number_of_halfedges() << std::endl;
        std::cout << "\t Number of facets    :\t" << surface_mesh.number_of_faces() << std::endl;
        PMP::stitch_borders(surface_mesh);
        std::cout << "Stitching done : " << std::endl;
        std::cout << "\t Number of vertices  :\t" << surface_mesh.number_of_vertices() << std::endl;
        std::cout << "\t Number of halfedges :\t" << surface_mesh.number_of_halfedges() << std::endl;
        std::cout << "\t Number of facets    :\t" << surface_mesh.number_of_faces() << std::endl;
    }

    int num_selfintersections() {
        bool intersecting = PMP::does_self_intersect(
               surface_mesh,
               PMP::parameters::vertex_point_map(get(CGAL::vertex_point, surface_mesh)));
        std::cout
          << (intersecting ? "There are self-intersections." : "There is no self-intersection.")
          << std::endl;
        std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
        PMP::self_intersections(surface_mesh, std::back_inserter(intersected_tris));
        std::cout << intersected_tris.size() << " pairs of triangles intersect." << std::endl;
    
        std::vector<std::pair<face_descriptor, face_descriptor>>::iterator it;
        for(it = intersected_tris.begin(); it != intersected_tris.end(); it++) {
            std::cout << it->first << std::endl;
            std::cout << it->second << std::endl;
            std::cout << std::endl;
            /* std::cout << it[0] << std::endl; */
        }
    }

    void remesh(double target_edge_length, int nb_iter, bool protect_edges) {
        std::vector<edge_descriptor> border;
        PMP::border_halfedges(
                faces(surface_mesh),
                surface_mesh,
                boost::make_function_output_iterator(halfedge2edge(surface_mesh, border)));
        PMP::split_long_edges(border, target_edge_length, surface_mesh);
    
        std::cout << "done." << std::endl;
        std::cout << "Start remeshing" << std::endl << " (" << num_faces(surface_mesh) << " faces)..." << std::endl;
      
        PMP::isotropic_remeshing(
                faces(surface_mesh),
                target_edge_length,
                surface_mesh,
                PMP::parameters::number_of_iterations(nb_iter)
                .protect_constraints(protect_edges));
    }

    void fill_holes() {
        unsigned int nb_holes = 0;
        BOOST_FOREACH(halfedge_descriptor h, halfedges(surface_mesh)) {
            if(is_border(h, surface_mesh)) {
                std::vector<face_descriptor> patch_facets;
                std::vector<vertex_descriptor> patch_vertices;

                bool success = CGAL::cpp11::get<0>(
                        PMP::triangulate_refine_and_fair_hole(
                            surface_mesh,
                            h,
                            std::back_inserter(patch_facets),
                            std::back_inserter(patch_vertices),
                            PMP::parameters::vertex_point_map(get(CGAL::vertex_point, surface_mesh)).geom_traits(K()))); 

                std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
                std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
                std::cout << "  Is fairing successful: " << success << std::endl;
                nb_holes++;
            }
        }
      
        std::cout << std::endl;
        std::cout << nb_holes << " holes have been filled" << std::endl;
    }
};
