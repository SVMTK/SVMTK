#include "CGALSurface.h"
#include "read_polygons_STL.h"
#include "surface_mesher.h"



// needed 

// TODO: Clean up PMP  header files
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/squared_distance_3.h> 

#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Mesh_3/Detect_features_in_polyhedra.h>
#include <CGAL/Poisson_reconstruction_function.h>

// CGAL IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <assert.h>
#include <boost/foreach.hpp>
#include <iterator>

#include <CGAL/IO/OFF_reader.h> // needed for load surface  

template< typename CGALSurface> // TODO: Conisder move to seperate header
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

template< typename Mesh> // seperate because -> typedef mesh in CGALSurface
bool load_surface(std::string filename, Mesh& mesh)
{
    std::ifstream input(filename);
    std::string file(filename);
    std::string extension = file.substr(file.find_last_of(".")+1);//TODO: FIX boost linking problem
    if (!input) {
        std::cerr << "Cannot open file " << std::endl;
        return false;
    }

    typedef typename Mesh::Point Point_3;
    std::vector<Point_3> points;
    std::vector< std::vector<std::size_t> > polygons;

    if ( extension=="off") {

       //if (input >> mesh) // NO ERROR, gives error { EROR
       //     return true;
       //}

       std::cout<< "reading off" << std::endl;
       if (!CGAL::read_OFF(input, points, polygons)) {
           std::cerr << "Error parsing the OFF file " << std::endl;
           return false;
       }
       std::cout<< "finished" << std::endl;
    }
    else if ( extension=="stl") {
       if (!read_polygons_STL(input, points, polygons)) { // TODO: the connection causes errors 
           std::cerr << "Error parsing the STL file " << std::endl;
           return false;
       }
    }
    else {
         std::cerr << "Error unkown file" << std::endl;
         return false;
    }
    std::cout<< "Step 1" << std::endl;
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
    std::cout<< "Step 2" << std::endl;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons,mesh);

    if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh))) {
       std::cout<< "reverse_face_orientation"<< std::endl;

       CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
    }
    return true;
}



//--------------------------------------------
//      CGALSURFACE FUNCTIONS
//--------------------------------------------
CGALSurface::Mesh& CGALSurface::get_mesh() {
    return mesh;
}

CGALSurface::CGALSurface(std::string filename) {
    load_surface(filename, mesh);
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

     // CGAL::Polygon_mesh_processing::detect_sharp_edges(mesh, 90);  included in cgal 4.13
     CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);

     CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(mesh),
                              target_edge_length,
                              mesh,
                              CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
                              .protect_constraints(protect_border));

}

int CGALSurface::fill_holes() {
    unsigned int nb_holes = 0;
    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh))
    {
      if(is_border(h, mesh))
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
void CGALSurface::insert_surface(CGALSurface& surface) // redunded 
{
  mesh+=surface.get_mesh();
}


CGALSurface::vertex_vector CGALSurface::points_inside(CGALSurface& other)
{
  CGALSurface::vertex_vector result;
  CGALSurface::Inside inside_poly2(other.get_mesh());

  for ( vertex_descriptor v_it : mesh.vertices())
  {
      CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
      if (res == CGAL::ON_BOUNDED_SIDE or res == CGAL::ON_BOUNDARY){result.push_back(v_it);}
  }
  return result;
}

CGALSurface::vertex_vector CGALSurface::points_outside(CGALSurface& other)
{
  CGALSurface::vertex_vector result;
  CGALSurface::Inside inside_poly2(other.get_mesh());

  for ( vertex_descriptor v_it : mesh.vertices())
  {
      CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
      if (res == CGAL::ON_UNBOUNDED_SIDE){result.push_back(v_it);}
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
  std::cout<< "here" << std::endl;
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
void CGALSurface::adjust_boundary(const double c)
{
    Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    CGALSurface::adjusting_boundary_region(vb, ve,c);
    //CGALSurface::adjusting_boundary_region(mesh.vertices().begin(), mesh.vertices.end(),c);
}
void CGALSurface::smooth_laplacian(const double c)
{
    Mesh::Vertex_range::iterator  vb = mesh.vertices().begin(), ve=mesh.vertices().end();
    CGALSurface::smooth_laplacian_region(vb, ve,c);
}

void CGALSurface::make_cube( double x0, double y0, double  z0,  double x1, double y1, double z1) // TODO :allow tilted cubes
{

        assert(x0!=x1);
        assert(y0!=y1);
        assert(z0!=z1);

        Mesh m;
        vertex_descriptor v0 = m.add_vertex(Point_3(x0,y0,z0));
        vertex_descriptor v1 = m.add_vertex(Point_3(x0,y0,z1));
        vertex_descriptor v2 = m.add_vertex(Point_3(x1,y0,z0));

        vertex_descriptor v3 = m.add_vertex(Point_3(x1,y0,z1));
        vertex_descriptor v4 = m.add_vertex(Point_3(x1,y1,z0));
        vertex_descriptor v5 = m.add_vertex(Point_3(x1,y1,z1));

        vertex_descriptor v6 = m.add_vertex(Point_3(x0,y1,z0));
        vertex_descriptor v7 = m.add_vertex(Point_3(x0,y1,z1));
  
        // Side 1
        m.add_face(v0, v2, v1);
        m.add_face(v3, v1, v2);

        // Side 2
        m.add_face(v3 ,v2, v5);
        m.add_face(v4 ,v5,v2);

        // Side 3 
        m.add_face(v4, v6, v5);
        m.add_face(v7, v5, v6);

        // Side 4
        m.add_face(v0, v1, v6);
        m.add_face(v7, v6, v1);

        // Side 5 
        m.add_face(v3, v5, v1);
        m.add_face(v7, v1, v5);

        // Side 6 
        m.add_face(v0, v6, v2);
        m.add_face(v4, v2, v6);

        this->mesh = m;

}

void CGALSurface::make_cylinder( double x0, double y0, double  z0,  double x1, double y1, double z1, double radius,  int number_of_segments) // TODO :allow tilted cubes
{
      
        CGALSurface::make_cone( x0, y0,z0, x1,y1, z1, radius , radius,  number_of_segments ) ;


}


void CGALSurface::make_cone( double x0, double y0, double  z0,  double x1, double y1, double z1, double r0, double r1,  int number_of_segments) // TODO :allow tilted cubes
{
        //TODO: Handle case for r1=0 or r0=0
        Mesh m;

    
        Index v0 = m.add_vertex(Point_3(x0,y0,z0));
        Index v1 = m.add_vertex(Point_3(x1,y1,z1));
        
        Index vb;
        Index vt;
        //normalizing vectors ...

        double n0 = x1-x0;
        double n1 = y1-y0;
        double n2 = z1-z0;

        double l1 = std::sqrt(n0*n0+n1*n1+n2*n2);

        Vector_3 normal( n0/l1 , n1/l1, n2/l1 );

        double l2 = std::sqrt( (n1-n2)*(n1-n2) + (n2-n0)*(n2-n0) + ( n0 -n1)*(n0-n1) );

        Vector_3 t1( (n1-n2)/l2 , (n2-n0)/l2, (n0-n1)/l2 ); 

        Vector_3 t2 = CGAL::cross_product(normal,t1); // normalized ?? 

        double c  = 360/number_of_segments;
     
        std::cout << normal << " " << t1 << " " << t2 << std::endl;
        for(int i = 0; i < number_of_segments; ++i) 
        {
               

           Point_3 pb= Point_3(x0,y0,z0)  + t1*r0*std::cos(c*i*CGAL_PI/180) + t2*r0*std::sin(c*i*CGAL_PI/180); 

           vb = m.add_vertex(pb);

           Point_3 pt = Point_3(x1,y1,z1)  + t1*r1*std::cos(c*i*CGAL_PI/180) + t2*r1*std::sin(c*i*CGAL_PI/180); 

           vt = m.add_vertex(pt);

           if ( i!=0)
           {  

             m.add_face(v0, vb,Index(vb-2));
             m.add_face(v1, Index(vt-2), vt);

             m.add_face(Index(vt-2), Index(vb-2), vt);
             m.add_face(vb, vt, Index(vb-2));

           }

        }
        m.add_face(Index(0), Index(2), vb);
        m.add_face(Index(1), vt, Index(3));

        m.add_face(vt, vb, Index(3));
        m.add_face(Index(2), Index(3), vb);



        this->mesh = m;
}

struct sphere_wrapper{
        
     public:
        static double radius;
        static double x0;
        static double y0;
        static double z0;
        
        static double function(double x, double y , double z) { return  (x-x0)*(x -x0) +  (y-y0)*(y -y0) + (z-z0)*(z -z0) -radius*radius; }

};

// TODO : Sphere uses a implicit function, and can be extended to other implict fucntion with smooth boundaries.  
double sphere_wrapper::radius = 0;
double sphere_wrapper::x0 = 0;
double sphere_wrapper::y0 = 0;
double sphere_wrapper::z0 = 0;

void CGALSurface::make_sphere( double x0, double y0, double  z0,double r0 ) 
{

  sphere_wrapper sphere;
   
  sphere.radius = r0;
  sphere.x0=x0;
  sphere.y0=y0;
  sphere.z0=z0;

  surface_mesher(mesh,sphere.function,x0,y0,z0,r0,30,r0/5,r0/5);   

}

void CGALSurface::split_edges(double  target_edge_length)
{
     CGAL::Polygon_mesh_processing::split_long_edges(edges(mesh), target_edge_length,mesh);
}


template<typename Implicit_function>
CGALSurface::CGALSurface(Implicit_function implicit_function,
             double bounding_sphere_radius,
             double angular_bound,
             double radius_bound,
             double distance_bound)
{
     surface_mesher(mesh,implicit_function,bounding_sphere_radius,angular_bound,radius_bound, distance_bound);
}

CGALSurface::CGALSurface()
{
   // constructor
}
void CGALSurface::save(const char* outpath)
{
  // TODO: save surface mesh as VTK.
    std::ofstream out(outpath);
    out << mesh;
    out.close();
}

void CGALSurface::preprocess(double target_edge_length,int nb_iter)
{
    // TODO: Implement deafult processing stream fix surface mesh

  
    //CGALSurface::triangulate_faces(mesh);
    //CGALSurface::isotropic_remeshing(mesh,target_edge_length,nb_iter,false);

}
void CGALSurface::fair(CGALSurface::vertex_vector vector)
{
     CGAL::Polygon_mesh_processing::fair(mesh,vector);
}





void CGALSurface::insert_points(std::vector<Point_3>& points)  // template iterator ?? ADD POINT ineffective, but could be implemented
{
     for( std::vector<Point_3>::iterator it=points.begin(); it!= points.end(); ++it){mesh.add_vertex(*it);}
}

void CGALSurface::poisson_reconstruction() 
{

    // TODO: Implement
    // clear mesh 
    //poisson_reconstruction();// all stored points , i.e. mesh. vertices_begin?? unpack to points
}




void CGALSurface::fix_close_junctures(double c)
{
   // TODO :  Clean up typedefs
   int count=0;
   std::cout<< "Step :" << ++count << std::endl;
        typedef boost::graph_traits<Mesh>::vertex_descriptor                     Point;
        typedef boost::graph_traits<Mesh>::vertices_size_type                    size_type;
        typedef boost::property_map<Mesh,CGAL::vertex_point_t>::type             Vertex_point_pmap;
        typedef CGAL::Search_traits_3<Kernel>                                    Traits_base;
        typedef CGAL::Search_traits_adapter<Point,Vertex_point_pmap,Traits_base> Traits;
        typedef CGAL::Orthogonal_k_neighbor_search<Traits>                       K_neighbor_search;
        typedef K_neighbor_search::Tree                                          Tree;
        typedef Tree::Splitter                                                   Splitter;
        typedef K_neighbor_search::Distance                                      Distance;

 
   CGALSurface::vertex_vector results;
   Vertex_point_pmap vppmap = get(CGAL::vertex_point,mesh);
   std::cout<< "Step :" << ++count << std::endl;
   Tree tree(vertices(mesh).begin(),
            vertices(mesh).end(),
            Splitter(),
            Traits(vppmap)
   );
   std::cout<< "Step :" << ++count << std::endl;
   Distance tr_dist(vppmap);

   // For each point in surface mesh -> find closest point if it is not a neighbor vertex -> add to vector
   FT distance;
   FT edgeL;
   Point_3 closest;
   bool flag;
   for (vertex_descriptor v_it : mesh.vertices())
   {
        flag = true;

        K_neighbor_search search(tree, mesh.point(v_it), 2,0,true,tr_dist); // DISTANCE CAN BE INCLUDED HEREmesh.point(v_it)
        
        closest = mesh.point((search.begin()+1)->first); //?? closest point is the query point, therefore chose second point 
        

        Point_3 current = mesh.point(v_it);
        distance = CGAL::squared_distance(current, closest );
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(v_it),mesh), done(vbegin);
        //ditr_dist.inverse_of_transformed_distance(search.begin()->second);
        do
        {
  
           edgeL = CGAL::squared_distance(current, mesh.point(*vbegin) ); 
           std::cout<< distance << " " << edgeL << std::endl;
           if ( distance >= edgeL ) 
           {
              flag=false;
            break;
           }
         *vbegin++;
        }while(vbegin!=done);

        if (flag){results.push_back(v_it);}

   }

 

   std::cout<< results.size() << std::endl;
   
   adjusting_boundary_region(results.begin(),results.end(),c);


   
}


CGALSurface::CGALSurface(const std::string f) {
    utils::read_off(mesh, f);
}


CGALSurface::CGALSurface(Polyhedron &polyhedron) {
    CGAL::copy_face_graph(polyhedron, mesh);
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
    Mesh::Vertex_range::iterator vb = mesh.vertices().begin(), ve = mesh.vertices().end();
    CGALSurface::adjusting_boundary_region(vb, ve, c);
}


void CGALSurface::smooth_laplacian(const double c) {
    Mesh::Vertex_range::iterator vb = mesh.vertices().begin(), ve = mesh.vertices().end();
    CGALSurface::smooth_laplacian_region(vb, ve, c);
}


void CGALSurface::smooth_taubin(const size_t nb_iter) {
    for (size_t i = 0; i < nb_iter; ++i) {
        this->smooth_laplacian(0.8);
        this->smooth_laplacian(-0.805);
    }
}


template<typename InputIterator>
void CGALSurface::adjusting_boundary_region(
        InputIterator begin,
        InputIterator end,
        const double c) {
    std::vector<std::pair<vertex_descriptor, Point_3> > smoothed; //rename
    for ( ; begin != end; ++begin) {
        Vector_3 delta = CGAL::Polygon_mesh_processing::compute_vertex_normal(*begin,mesh);
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
        CGAL::Vertex_around_target_circulator<Mesh> vbegin(mesh.halfedge(*begin), mesh), done(vbegin);
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


int CGALSurface::num_self_intersections() {
    std::vector< std::pair<face_descriptor, face_descriptor> > intersected_tris;
    CGAL::Polygon_mesh_processing::self_intersections(mesh, std::back_inserter(intersected_tris));
    return intersected_tris.size();  // Could actually return the triangles themselves
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


std::vector<vertex_descriptor> CGALSurface::points_inside(CGALSurface &other) {
    std::vector<vertex_descriptor> result;
    inside inside_poly2(other.get_mesh());
    for (const auto &v_it: mesh.vertices()) {
        const CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
        if (res == CGAL::ON_BOUNDED_SIDE or res == CGAL::ON_BOUNDARY) {
            result.push_back(v_it);
        }
    }
    return result;
}


std::vector<vertex_descriptor> CGALSurface::points_outside(CGALSurface& other) {
    std::vector<vertex_descriptor> result;
    inside inside_poly2(other.get_mesh());

    for (const auto &v_it: mesh.vertices()) {
        const CGAL::Bounded_side res = inside_poly2(mesh.point(v_it));
        if (res == CGAL::ON_UNBOUNDED_SIDE) {
              result.push_back(v_it);
        }
    }
    return result;
}


void CGALSurface::fair(std::vector<vertex_descriptor> &vector) {
    CGAL::Polygon_mesh_processing::fair(mesh, vector);
}


void CGALSurface::surface_intersection(CGALSurface &other) {
    CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh, other.get_mesh(), mesh);
}


void CGALSurface::surface_union(CGALSurface &other) {       // Probably bad to use union name
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(
            mesh,
            other.get_mesh(),
            mesh);
}


void CGALSurface::surface_difference(CGALSurface &other) {
    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(mesh, other.get_mesh(), mesh);
}


void CGALSurface::reconstruct_surface(
        const double sm_angle,
        const double sm_radius,
        const double sm_distance) {
    RPolyhedron input_polyhedron;
    this->get_polyhedron(input_polyhedron);
    std::deque<RPwn> points;

    // TODO: Can I do some sort of casting to avoid copying the whole mesh above?
    for (boost::graph_traits<RPolyhedron>::vertex_descriptor vd: vertices(input_polyhedron)) {
        const RPoint p = vd->point();
        const RVector n = CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, input_polyhedron);
        points.push_back(std::make_pair(p, n));
    }

    RPolyhedron output_mesh;

    double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag> (points, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>()));

    CGAL::poisson_surface_reconstruction_delaunay(
            points.begin(), points.end(),
            CGAL::First_of_pair_property_map<RPwn>(),
            CGAL::Second_of_pair_property_map<RPwn>(),
            output_mesh, average_spacing);

    CGAL::copy_face_graph(output_mesh, this->mesh);
}


int CGALSurface::num_faces() const {
    return mesh.number_of_faces();
}


int CGALSurface::num_edges() const {
    return mesh.number_of_edges();
}


int CGALSurface::num_vertices() const {
    return mesh.number_of_vertices();
}













