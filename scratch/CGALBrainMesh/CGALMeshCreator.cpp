#include "CGALMeshCreator.h"
#include "remove_isolated_vertices.h"


CGALMeshCreator::CGALMeshCreator(std::vector<CGALSurface> surfaces, CGAL::Bbox_3 bbox_3, AbstractMap& map)
{
    Function_vector v;
    Polyhedron polyhedron;
    for(std::vector<CGALSurface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       sit->get_polyhedron(polyhedron);
       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       v.push_back(polyhedral_domain);
    }
    Function_wrapper wrapper(v,map);
    Mesh_domain domain(wrapper,bbox_3);
    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,bbox_3));
    default_parameters();
}

CGALMeshCreator::CGALMeshCreator( std::vector<CGALSurface> surfaces , AbstractMap& map )
{
    Function_vector v;
    Polyhedron polyhedron;
    for(std::vector<CGALSurface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       sit->get_polyhedron(polyhedron);
       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       v.push_back(polyhedral_domain);
    }
    Function_wrapper wrapper(v,map);
    Mesh_domain domain(wrapper,wrapper.bbox());

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,wrapper.bbox()));
    default_parameters();
}

CGALMeshCreator::CGALMeshCreator( CGALSurface& surface, CGAL::Bbox_3 bbox_3)
{
    Function_vector v;
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);
    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
    v.push_back(polyhedral_domain);
    Function_wrapper wrapper(v);
    Mesh_domain domain(wrapper,bbox_3);
    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,bbox_3));
    default_parameters();
}

CGALMeshCreator::CGALMeshCreator(CGALSurface &surface)
{
    Function_vector v;
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);
    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
    v.push_back(polyhedral_domain);
    Function_wrapper wrapper(v);
    Mesh_domain domain(wrapper,wrapper.bbox());
    default_parameters();

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,wrapper.bbox()));
}



template < typename InputIterator>
void CGALMeshCreator::add_polylines(InputIterator begin, InputIterator end)
{
    domain_ptr->add_features(begin, end);
}
void CGALMeshCreator::add_polylines(const char* filename)
{
    typedef std::vector<Point_3> Polyline_3;

    //std::vector<Polyline_3>& polylines;
    //domain_ptr->add_features(polylines.begin(), polylines.end());


}


void CGALMeshCreator::lipschitz_size_field(int subdomain_id, int k,double min_size,double max_size)
{
    if (lip_sizing_ptr)
    {
       lip_sizing_ptr->add_parameters_for_subdomain(subdomain_id,k,min_size,max_size);
    }
    else
    {
       lip_sizing_ptr = std::unique_ptr<Lip_sizing> (new Lip_sizing(*domain_ptr.get()));
       lip_sizing_ptr.get()->add_parameters_for_subdomain(subdomain_id,k,min_size,max_size);

    }
}

void CGALMeshCreator::set_parameters(Parameters new_parameters)
{
    for (Parameters::iterator pit= new_parameters.begin(); pit!=new_parameters.end(); ++pit )
    {
        parameters[pit->first] = static_cast<double>(pit->second);
        std::cout << pit->first <<" " << pit->second << std::endl;
    }
}
void CGALMeshCreator::set_parameter(std::string key , double value )
{
   parameters[key] = value;
}




void CGALMeshCreator::create_mesh()
{
    std::cout << "begin_meshing" << std::endl;
    Mesh_criteria criteria(CGAL::parameters::edge_size  =parameters["edge_size"],
                           CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );


    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria);

    std::cout << "done_meshing" << std::endl;


//    if(parameters["odt_optimize"] ) do it manually?
//    {
//    }
//    if(parameters["perturb_optimize"])
//    {
//    }
//    if(parameters["lloyd_optimize"])
//    {
//    }
//    if( parameters["exude_optimize"] )
//    {
//    }
    //remove_isolated_vertices(c3t3);


}
void CGALMeshCreator::create_mesh(int inital_points)
{
    Mesh_criteria criteria(CGAL::parameters::edge_size  =parameters["edge_size"],
                           CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );

    CGAL::internal::Mesh_3::init_c3t3(c3t3, *domain_ptr.get(), criteria, inital_points);

    refine_mesh_3(c3t3, *domain_ptr.get(), criteria,CGAL::parameters::no_reset_c3t3());


//    if(parameters["odt_optimize"] )
//    {

//    }

//    if(parameters["perturb_optimize"])
//    {

//    }
//    if(parameters["lloyd_optimize"])
//    {

//    }
//    if( parameters["exude_optimize"] )
//    {

//    }
    remove_isolated_vertices(c3t3);

}

void CGALMeshCreator::save_mesh(const char* OutPath)
{
    std::ofstream  medit_file(OutPath);
    c3t3.output_to_medit(medit_file);
    medit_file.close();
}
