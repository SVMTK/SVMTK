#ifndef __CGAL_MESH_CREATOR_H
#define __CGAL_MESH_CREATOR_H


// LOCAL
#include "CGALSurface.h"
#include "SubdomainMap.h"
#include "Polyhedral_vector_to_labeled_function_wrapper.h"

// STD
#include <list>
#include <fstream>
#include <memory>

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_3/experimental/Lipschitz_sizing_polyhedron.h>

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/remove_far_points_in_mesh_3.h>
#include <CGAL/config.h>
#include <CGAL/assertions.h>

#include <CGAL/refine_mesh_3.h>

#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h> // FIX

#include <CGAL/IO/File_medit.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Mesh_3/Detect_polylines_in_polyhedra.h>
#include <CGAL/Mesh_3/polylines_to_protect.h>


template<typename C3T3>
int remove_isolated_vertices(C3T3 &c3t3) {
    typedef typename C3T3::Triangulation Tr;
    /* typedef typename C3T3::Cells_in_complex_iterator Cell_iterator; */

    /* typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator; */
    typedef typename Tr::Vertex_handle Vertex_handle;

    std::map<Vertex_handle, bool> vertex_map;
    for( auto vit = c3t3.triangulation().finite_vertices_begin();
            vit != c3t3.triangulation().finite_vertices_end();
            ++vit)
    {
        vertex_map[vit] = false;
    }

    for(auto cit = c3t3.cells_in_complex_begin(); cit != c3t3.cells_in_complex_end(); ++cit)
    {
        for (std::size_t i = 0; i < 4; ++i)
            vertex_map[cit->vertex(i)] = true;
    }

    const auto before = c3t3.triangulation().number_of_vertices();
    for (auto &it: vertex_map)
    {
        if (!it.second)
            c3t3.triangulation().remove(it.first);
    }
    const auto after = c3t3.triangulation().number_of_vertices();
    std::cout << "Number of vertices removed: " << before - after << std::endl;
    return before - after;
}


class CGALSurface;

class CGALMeshCreator {
    public:

        //---------------------------------------------------------------
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef K::Point_3 Point_3;

        typedef CGAL::Mesh_polyhedron_3< K >::type Polyhedron;

        typedef CGAL::Polyhedral_mesh_domain_with_features_3< K, Polyhedron > Polyhedral_mesh_domain_3;
        typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper< Polyhedral_mesh_domain_3, K > Function_wrapper;

        typedef Function_wrapper::Function_vector Function_vector;
        typedef CGAL::Labeled_mesh_domain_3< K> Labeled_Mesh_Domain;
        typedef CGAL::Mesh_domain_with_polyline_features_3< Labeled_Mesh_Domain > Mesh_domain;

        typedef CGAL::Mesh_3::Lipschitz_sizing< K, Mesh_domain > Lip_sizing;
        typedef CGAL::Mesh_triangulation_3< Mesh_domain >::type Tr;

        typedef Mesh_domain::Curve_segment_index Curve_index;
        typedef Mesh_domain::Corner_index Corner_index;

        typedef CGAL::Mesh_complex_3_in_triangulation_3< Tr, Corner_index, Curve_index > C3t3;

        typedef CGAL::Mesh_criteria_3< Tr > Mesh_criteria;
        typedef Tr::Vertex_handle Vertex_handle;
        /* typedef Tr::Finite_vertices_iterator Finite_vertices_iterator; */
        /* typedef Tr::Locate_type Locate_type; */

        typedef C3t3::Subdomain_index Subdomain_index;
        /* typedef C3t3::Cells_in_complex_iterator Cell_iterator; */

        typedef std::vector< Point_3 > Polyline_3;
        typedef std::vector< Polyline_3 > Polylines;
        typedef std::map< std::string, double > Parameters;


        struct Minimum_sphere {
            typedef CGAL::Min_sphere_of_spheres_d_traits_3<K, K::FT> MinSphereTraits;
            typedef CGAL::Min_sphere_of_spheres_d<MinSphereTraits> Min_sphere;
            typedef MinSphereTraits::Sphere Sphere;

            template<typename MeshPolyhedron_3>
            void add_polyhedron(const MeshPolyhedron_3 &polyhedron) {
                for (auto it: vertices(polyhedron)) {
                    S.push_back(Sphere(it->point(), 0.0));
                }
            }

            double get_bounding_sphere_radius() {
                Min_sphere ms(S.begin(), S.end());
                return CGAL::to_double(ms.radius());
            }

            private:
                std::vector<Sphere> S;
        };


        CGALMeshCreator(CGALSurface& surface);
        CGALMeshCreator(std::vector<CGALSurface> &surfaces);
        CGALMeshCreator(std::vector<CGALSurface> &surfaces, AbstractMap &map);

        ~CGALMeshCreator() {}


        /* // TODO: Implement volume ++ */ 

        void set_parameters(const Parameters &new_parameters);
        void set_parameter(const std::string key, const double value) { parameters[key] = value; };
        void default_parameters() {
            parameters["mesh_resolution"] = 64.0;
            parameters["facet_angle"] = 25.0;
            parameters["facet_size"] = 0.1;
            parameters["facet_distance"] = 0.1;
            parameters["cell_radius_edge_ratio"] = 3.0;
            parameters["cell_size"] = 0.1;
            parameters["edge_size"] = 0.1;
        }

        void default_creating_mesh();       // Remove or set as default?
        void create_mesh();
        void create_mesh(const double mesh_resolution);

        void refine_mesh(const double mesh_resolution );
        void refine_mesh();

        void save(const std::string outpath);

        Polylines& get_features() {return features; }
        void add_feature(Polyline_3 polyline) {features.push_back(polyline); }


        template<typename InputIterator>
        void set_features(InputIterator begin, InputIterator end){domain_ptr->add_features(begin, end); }
        void set_features(Polylines& polylines) {set_features(polylines.begin(), polylines.end() ); }
        void set_features();

        Polylines& get_borders() {return borders; }
        void set_borders() {domain_ptr.get()->add_features(get_borders().begin(), get_borders() .end());}
        void add_borders(Polyline_3 polyline) { borders.push_back(polyline); }

        void add_sharp_border_edges(Polyhedron& polyhedron);
        void add_sharp_border_edges(CGALSurface& surface);

        void reset_borders(){ borders.clear(); }

        std::shared_ptr<CGALSurface> get_boundary();

        void lloyd(const double time_limit = 0, const int max_iteration_number = 0,
                const double convergence = 0.02, const double freeze_bound = 0.01,
                const bool do_freeze = true);

        void odt(const double time_limit= 0, const int max_iteration_number = 0,
                const double convergence = 0.02, const double freeze_bound = 0.01,
                const bool do_freeze = true);

        void excude(const double tl = 0, const double sb = 0) { CGAL::exude_mesh_3(c3t3, sb, tl); }

        void perturb(const double tl = 0, const double sb = 0) { CGAL::perturb_mesh_3 ( c3t3, *domain_ptr.get(), tl, sb); }

        void lipschitz_size_field(const int subdomain_id, const int k, const double min_size, const double max_size);

        void label_boundary_cells(const int btag, const int ntag);
        void remove_label_cells(const int tag);

        int number_of_cells() { return c3t3.number_of_cells(); }

    private:
        std::unique_ptr<Mesh_domain> domain_ptr;
        std::unique_ptr<Lip_sizing> lip_sizing_ptr;
        Minimum_sphere min_sphere;

        Polylines borders;
        Polylines features;
        Parameters parameters;
        C3t3 c3t3;
};


CGALMeshCreator::CGALMeshCreator(std::vector<CGALSurface> &surfaces)
{
    /**********************************
    *  Constructor for CGALSurfaces
    **********************************/
    Function_vector v;
    Polyhedron polyhedron;

    for(auto &sit: surfaces)
    {
       sit.get_polyhedron(polyhedron);
       min_sphere.add_polyhedron(polyhedron);
       v.push_back(new Polyhedral_mesh_domain_3(polyhedron));
    }

    Function_wrapper wrapper(v);
    domain_ptr = std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,wrapper.bbox()));
    default_parameters();
}


CGALMeshCreator::CGALMeshCreator(std::vector<CGALSurface> &surfaces, AbstractMap &amap)
{
    /**********************************
    *  Constructor for CGALSurfaces
    **********************************/
    Function_vector v;
    Polyhedron polyhedron;

    for(auto &sit: surfaces)
    {
        sit.get_polyhedron(polyhedron);
        min_sphere.add_polyhedron(polyhedron);
        v.push_back(new Polyhedral_mesh_domain_3(polyhedron));
    }

    Function_wrapper wrapper(v, amap);
    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));
    default_parameters();
}


CGALMeshCreator::CGALMeshCreator(CGALSurface &surface)
{
    /**********************************
    *  Constructor for CGALSurfaces
    **********************************/
    Function_vector v;
    Polyhedron polyhedron;

    surface.get_polyhedron(polyhedron);
    min_sphere.add_polyhedron(polyhedron);

    v.push_back(new Polyhedral_mesh_domain_3(polyhedron));
    Function_wrapper wrapper(v);
    domain_ptr = std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));
    default_parameters();
}


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------


void CGALMeshCreator::set_parameters(const Parameters &new_parameters) {
    for (auto pit: new_parameters) {
        parameters[pit.first] = static_cast<double>(pit.second);
    }
}


void CGALMeshCreator::lipschitz_size_field(
        const int subdomain_id, const int k, const double min_size, const double max_size)
{
    // lipschitz_size_field  
    if (lip_sizing_ptr) {
       lip_sizing_ptr->add_parameters_for_subdomain(subdomain_id, k, min_size, max_size);
    }
    else {
        lip_sizing_ptr = std::unique_ptr<Lip_sizing> (new Lip_sizing(*domain_ptr.get()));
        lip_sizing_ptr.get()->add_parameters_for_subdomain(subdomain_id, k, min_size, max_size);
    }
}


void CGALMeshCreator::create_mesh()
{
    Mesh_criteria criteria(CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );

    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);
    remove_isolated_vertices(c3t3);
}


void CGALMeshCreator::create_mesh(const double mesh_resolution)
{
    const auto r = min_sphere.get_bounding_sphere_radius();
    const auto cell_size = r/mesh_resolution*2.0;

    Mesh_criteria criteria(CGAL::parameters::edge_size = cell_size,
                                       CGAL::parameters::facet_angle = 30.0,
                                       CGAL::parameters::facet_size = cell_size,
                                       CGAL::parameters::facet_distance = cell_size/10.0, 
                                       CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                       CGAL::parameters::cell_size = cell_size);

    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);
    remove_isolated_vertices(c3t3);
}


void CGALMeshCreator::default_creating_mesh()
{
    // TODO: Is this necessary
    Mesh_criteria criteria(CGAL::parameters::facet_angle = 25.0,
                           CGAL::parameters::edge_size = 0.025,
                           CGAL::parameters::facet_size = 0.05,
                           CGAL::parameters::facet_distance = 0.005,
                           CGAL::parameters::cell_radius_edge_ratio = 3,
                           CGAL::parameters::cell_size = 0.05);

    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);
    remove_isolated_vertices(c3t3);
}


void CGALMeshCreator::save(const std::string outpath)
{
    std::ofstream medit_file(outpath);
    c3t3.output_to_medit(medit_file);
    medit_file.close();
}


//----------------------   Have overloaded functions or just one?
void CGALMeshCreator::add_sharp_border_edges(Polyhedron& polyhedron) {
    Polylines polylinput;

    auto eif = get(CGAL::edge_is_feature, polyhedron);
    CGAL::Polygon_mesh_processing::detect_sharp_edges(polyhedron, 80, eif); // -> threshold ?? 
    for(auto he = polyhedron.edges_begin(); he != polyhedron.edges_end(); ++he) {
        if(he->is_feature_edge()) {
            Polyline_3 polyline;
            polyline.push_back(he->vertex()->point());
            polyline.push_back(he->opposite()->vertex()->point());
            polylinput.push_back(polyline);
        }
    }
    polylines_to_protect(this->borders, polylinput.begin(), polylinput.end());
}


void CGALMeshCreator::add_sharp_border_edges(CGALSurface& surface)
{
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);
    Polylines polylinput;

    auto eif = get(CGAL::edge_is_feature, polyhedron);
    CGAL::Polygon_mesh_processing::detect_sharp_edges(polyhedron, 80, eif);
    for (auto he = polyhedron.edges_begin(); he != polyhedron.edges_end(); ++he)
    {
        if(he->is_feature_edge())
        {
            Polyline_3 polyline;
            polyline.push_back(he->vertex()->point());
            polyline.push_back(he->opposite()->vertex()->point());
            polylinput.push_back(polyline);
        }
    }
    polylines_to_protect(this->borders, polylinput.begin(), polylinput.end());
}


void CGALMeshCreator::set_features()
{
    Polylines polylines;
    polylines_to_protect(polylines, features.begin(), features.end());
    this->set_features(polylines.begin(), polylines.end());
}


void CGALMeshCreator::lloyd(const double tl, const int max_iter, const double conv,
        const double fb, const bool df)
{
    CGAL::lloyd_optimize_mesh_3(c3t3, *domain_ptr.get(), tl, max_iter, conv, fb, df);
}


void CGALMeshCreator::odt(const double tl, const int max_iter, const double conv,
        const double fb, const bool df)
{
    CGAL::odt_optimize_mesh_3(c3t3, *domain_ptr.get(), tl, max_iter, conv, fb, df);
}


void CGALMeshCreator::label_boundary_cells(const int btag, const int ntag)
{
    // workaround to mark boundary cells of for example lateral ventircles. this allows for easy marking of Facetfunction in FEniCS
    Subdomain_index subdomain_index(btag);
    Subdomain_index subdomain_index_bis(ntag);
    for(auto cit = c3t3.cells_in_complex_begin(subdomain_index); cit != c3t3.cells_in_complex_end(); ++cit)
    {
        for (std::size_t i = 0; i < 4; ++i)
        {
            if (c3t3.subdomain_index(cit->neighbor(i)) != subdomain_index)
                c3t3.set_subdomain_index(cit->neighbor(i), subdomain_index_bis);
        }
    }
}


void CGALMeshCreator::remove_label_cells(const int tag)
{
    Subdomain_index subdomain_index(tag);
    for(auto cit = c3t3.cells_in_complex_begin(subdomain_index); cit != c3t3.cells_in_complex_end(); ++cit) {
        c3t3.remove_from_complex(cit);
    }
}


void CGALMeshCreator::refine_mesh()
{
    // TODO: Use default parameters
    Mesh_criteria criteria(CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );

    refine_mesh_3(c3t3, *domain_ptr.get(), criteria,CGAL::parameters::no_reset_c3t3());

    while (remove_isolated_vertices(c3t3) > 0)
            c3t3.rescan_after_load_of_triangulation();
}


void CGALMeshCreator::refine_mesh(const double mesh_resolution)
{
    double r = min_sphere.get_bounding_sphere_radius();
    const double cell_size = r/mesh_resolution;

    Mesh_criteria criteria(CGAL::parameters::edge_size = cell_size,
                                       CGAL::parameters::facet_angle = 30.0,
                                       CGAL::parameters::facet_size = cell_size,
                                       CGAL::parameters::facet_distance = cell_size/10.0, 
                                       CGAL::parameters::cell_radius_edge_ratio = 3.0,
                                       CGAL::parameters::cell_size = cell_size);

   refine_mesh_3(c3t3, *domain_ptr.get(), criteria,CGAL::parameters::no_reset_c3t3());

   if ( remove_isolated_vertices(c3t3) > 0)
       c3t3.rescan_after_load_of_triangulation();
}


std::shared_ptr<CGALSurface> CGALMeshCreator::get_boundary()
{

   std::shared_ptr<CGALSurface> surf(new CGALSurface());
   if (c3t3.number_of_cells()> 0)
   {
        facets_in_complex_3_to_triangle_mesh(c3t3, surf->get_mesh());
   }
   else
   {
       std::cout << "No mesh generated, output is an empty surface file" << std::endl;
   }
   return surf;
}


#endif
