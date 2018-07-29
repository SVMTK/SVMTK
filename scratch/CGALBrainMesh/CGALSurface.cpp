#include "CGALSurface.h"
#include "read_polygons_STL.h"
//#include "surface_mesher.h"

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/IO/OFF_reader.h>

#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <boost/foreach.hpp>
#include <iterator>


template< typename CGALSurface>
void surface_overlapp(CGALSurface& surf1 , CGALSurface& surf2 )
{
    typedef typename CGALSurface::vertex_vector vertex_vector;

    vertex_vector surf1points;
    vertex_vector surf2points;

    surf1points = surf1.points_inside(surf2);
    surf2points = surf2.points_inside(surf1);
    double c=-0.1;

    while (  !surf1points.empty() and  !surf2points.empty() )
    {

            surf1.fair(surf1points);
            surf2.fair(surf2points);

            surf1.adjusting_boundary_region(surf1points.begin() ,surf1points.end(), c);
            surf2.adjusting_boundary_region(surf2points.begin() ,surf2points.end(), c);

            surf1points.clear();
            surf2points.clear();

            surf1points= surf1.points_inside(surf2);
            surf2points = surf2.points_inside(surf1);
    }

}

template< typename Mesh>
bool load_surface(char *filename, Mesh& mesh)
{
  std::ifstream input(filename);
  std::string file(filename);
  std::string extension = file.substr(file.find_last_of(".")+1);//TODO: FIX boost linking problem
  if (!input)
  {
    std::cerr << "Cannot open file " << std::endl;
    return false;
  }

  typedef typename Mesh::Point Point_3;
  std::vector<Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;

  if ( extension=="off")
  {
     if (!CGAL::read_OFF(input, points, polygons))
     {
         std::cerr << "Error parsing the OFF file " << std::endl;
         return false;
     }
  }
  else if ( extension=="stl")
  {
     if (!read_polygons_STL(input, points, polygons)) // TODO: the connection causes errors
     {
         std::cerr << "Error parsing the STL file " << std::endl;
         return false;
     }
  }
  else
  {
     std::cerr << "Error unkown file" << std::endl;
     return false;
  }

  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,mesh);

  if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
  {
   std::cout<< "reverse_face_orientation"<< std::endl;

    CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
  }
  return true;
}

CGALSurface::Mesh& CGALSurface::get_mesh()
{
    return mesh;
}

CGALSurface::CGALSurface(char *filename)
{
    load_surface(filename,mesh);
}

void CGALSurface::operator^= (CGALSurface& other )
{
   CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh, other.get_mesh(),mesh) ;
}

void CGALSurface::operator-=(CGALSurface& other)
{
    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(mesh , other.get_mesh(),mesh);
}

void CGALSurface::operator+=( CGALSurface& other )
{
   CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh, other.get_mesh(),mesh);
}
bool CGALSurface::self_intersections()
{
    return CGAL::Polygon_mesh_processing::does_self_intersect(mesh,
                CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));
}

void CGALSurface::isotropic_remeshing(double target_edge_length, unsigned int nb_iter, bool protect_border)
{

     CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);

     CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh),
                              target_edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
                              .protect_constraints(protect_border));

}

int CGALSurface::fill_holes()
{
    unsigned int nb_holes = 0;
    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh))
    {
      if(is_border(h,mesh))
      {
        std::vector<face_descriptor>  patch_facets;
        std::vector<vertex_descriptor> patch_vertices;
        bool success = CGAL::cpp11::get<0>(
          CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
                    mesh,
                    h,
                    std::back_inserter(patch_facets),
                    std::back_inserter(patch_vertices),
          CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).
                    geom_traits(Kernel())) );

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
bool CGALSurface::triangulate_faces()
{
    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

    BOOST_FOREACH(face_descriptor fit, faces(mesh))
      if (next(next(halfedge(fit, mesh), mesh), mesh)
          !=   prev(halfedge(fit, mesh), mesh))
        std::cerr << "Error: non-triangular face left in mesh." << std::endl;
    return true;
}
void CGALSurface::stitch_borders()
{
    CGAL::Polygon_mesh_processing::stitch_borders(mesh);

}
void CGALSurface::insert_surface(CGALSurface& surface)
{
  // TODO:
}


CGALSurface::vertex_vector CGALSurface::points_inside(CGALSurface& other)
{
  vertex_vector result;
  CGALSurface::Inside inside_poly2(other.get_mesh());

  for ( vertex_descriptor v_it : mesh.vertices())
  {
      std::cout << "ehi" << std::endl;
      CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
      if (res == CGAL::ON_BOUNDED_SIDE or res == CGAL::ON_BOUNDARY)
      {
            result.push_back(v_it);
      }
  }
  return result;
}

CGALSurface::vertex_vector CGALSurface::points_outside(CGALSurface& other)
{
  vertex_vector result;
  CGALSurface::Inside inside_poly2(other.get_mesh());

  for ( vertex_descriptor v_it : mesh.vertices())
  {
      CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
      if (res == CGAL::ON_UNBOUNDED_SIDE)
      {
            result.push_back(v_it);
      }
  }
  return result;
}



template<typename InputIterator >
void CGALSurface::adjusting_boundary_region(InputIterator begin , InputIterator  end, const double c)
{

  typedef typename CGAL::Vertex_around_target_circulator<Mesh> HV_const_circulator;
  std::vector<std::pair<vertex_descriptor, Point_3> > smoothed; //rename


  for ( ; begin != end; ++begin)
  {
      Vector_3 delta= CGAL::Polygon_mesh_processing::compute_vertex_normal(*begin,mesh);
      //delta = delta/std::sqrt(delta*delta);
      Point_3 p = mesh.point(*begin) + c*delta;
      smoothed.push_back(std::make_pair(*begin, p));
  }
  for (const std::pair<vertex_descriptor, Point_3> s : smoothed)
  {
      mesh.point(s.first) = s.second;
  }
}
template<typename InputIterator >
void CGALSurface::smooth_laplacian_region(InputIterator begin , InputIterator end ,const int c)
{
  typedef typename CGAL::Vertex_around_target_circulator<Mesh> HV_const_circulator;
  std::vector<std::pair<vertex_descriptor, Point_3> > smoothed;
  for ( ; begin != end; ++begin)
  {

      Point_3 current = mesh.point(*begin);
      Vector_3 delta=CGAL::NULL_VECTOR;
      CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(*begin),mesh), done(vbegin);
      do
      {
          delta += Vector_3(mesh.point(*vbegin) - current);
          *vbegin++;
      }while(vbegin!=done);

      //delta = delta/std::sqrt(delta*delta);
      Point_3 p = current + c*delta/mesh.degree(*begin);
      smoothed.push_back(std::make_pair(*begin, p));
  }
  for (const std::pair<vertex_descriptor, Point_3> s : smoothed)
  {
      mesh.point(s.first) = s.second;
  }


}
void CGALSurface::adjust_boundary(const double c) {
    Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    CGALSurface::adjusting_boundary_region(vb, ve,c);
    //CGALSurface::adjusting_boundary_region(mesh.vertices().begin(), mesh.vertices.end(),c);
}
void CGALSurface::smooth_laplacian(const double c)
{
    Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    CGALSurface::smooth_laplacian_region(vb, ve,c);
}


template< typename Polyhedron_3>
void CGALSurface::get_polyhedron(Polyhedron_3& polyhedron_3) {
    CGAL::copy_face_graph(mesh, polyhedron_3);
}


template<typename Implicit_function>
CGALSurface::CGALSurface(Implicit_function implicit_function,
             double bounding_sphere_radius,
             double angular_bound,
             double radius_bound,
             double distance_bound)
{
    // TODO: FIXME
     // surface_mesher(mesh,implicit_function,bounding_sphere_radius,angular_bound,radius_bound, distance_bound);
}
void CGALSurface::save(const char* outpath)
{
    std::ofstream out(outpath);
    out << mesh;
    out.close();
}

void CGALSurface::preprocess(double target_edge_length,int nb_iter)
{
    CGALSurface::triangulate_faces(mesh);
    CGALSurface::isotropic_remeshing(mesh,target_edge_length,nb_iter,false);

}
void CGALSurface::fair(CGALSurface::vertex_vector vector)
{
     CGAL::Polygon_mesh_processing::fair(mesh,vector);
}
