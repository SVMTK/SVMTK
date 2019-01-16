#include "CGALMeshCreator.h"
#include "remove_isolated_vertices.h"

// move from h other stuff
#include <CGAL/IO/File_medit.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// #include <CGAL/Mesh_3/Detect_features_in_polyhedra.h>

// Current problem simple geometries with sharp edges 90 degree 
// CLEAN : add egdes and insert edges
// insert_edges -> insert points, add_edges -> insert polylines

CGALMeshCreator::CGALMeshCreator(std::vector<CGALSurface> surfaces, CGAL::Bbox_3 bbox_3, AbstractMap& map)
{
    Function_vector v;
    Polyhedron polyhedron;
    for(std::vector<CGALSurface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       sit->get_polyhedron(polyhedron);
       // TODO: FIXME: FIX ME : detect features error 
       /* detect_features(polyhedron, 90); */
       /* add_sharp_edges(polyhedron); */

       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       v.push_back(polyhedral_domain);
    }
    Function_wrapper wrapper(v,map);
    Mesh_domain domain(wrapper,bbox_3);
    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,bbox_3));
    default_parameters();
}

CGALMeshCreator::CGALMeshCreator( std::vector<CGALSurface> surfaces )
{
    Function_vector v;
    Polyhedron polyhedron;
    for(std::vector<CGALSurface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       sit->split_edges(0.4);
       sit->get_polyhedron(polyhedron);
       // TODO: FIX ME
       /* detect_features(polyhedron, 90); */
       /* add_sharp_edges(polyhedron); */

       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       v.push_back(polyhedral_domain);
    }
    Function_wrapper wrapper(v);
    Mesh_domain domain(wrapper, wrapper.bbox());

    domain_ptr = std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,wrapper.bbox()));
    default_parameters();
}

CGALMeshCreator::CGALMeshCreator( std::vector<CGALSurface> surfaces , AbstractMap& map )
{
    Function_vector v;
    Polyhedron polyhedron;
    for(std::vector<CGALSurface>::iterator sit= surfaces.begin() ;sit!= surfaces.end();sit++)
    {
       sit->get_polyhedron(polyhedron);
       // TODO: FIX ME
       /* detect_features(polyhedron, 90); */
       /* add_sharp_edges(polyhedron); */

       Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
       v.push_back(polyhedral_domain);
    }


    Function_wrapper wrapper(v,map);
    Mesh_domain domain(wrapper,wrapper.bbox());

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper, wrapper.bbox()));
    default_parameters();
}


CGALMeshCreator::CGALMeshCreator( CGALSurface& surface, CGAL::Bbox_3 bbox_3)
{
    Function_vector v;
    Polyhedron polyhedron;
    surface.get_polyhedron(polyhedron);

    /* detect_features(polyhedron, 90); */
    /* add_sharp_edges(polyhedron); */

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


    // detect all  features in polyhedron : 
    /* detect_features(polyhedron, 90); */
    /* add_sharp_edges(polyhedron); */

    Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(polyhedron);
    v.push_back(polyhedral_domain);
    Function_wrapper wrapper(v);

    Mesh_domain domain(wrapper,wrapper.bbox());
   
    default_parameters();

    domain_ptr=std::unique_ptr<Mesh_domain> (new Mesh_domain(wrapper,wrapper.bbox()));

}


void CGALMeshCreator::lipschitz_size_field(int subdomain_id, int k,double min_size,double max_size)
{

    // lipschitz_size_field  
    // TODO : Explain 

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


    // TODO : Minimal  inital resticitons = None ? , to many restrictions causes segmentation fault



    //Mesh_criteria criteria(CGAL::parameters::edge_size  =parameters["edge_size"],
    //                       CGAL::parameters::facet_angle=parameters["facet_angle"],
    //                       CGAL::parameters::facet_size =parameters["facet_size"],
    //                       CGAL::parameters::facet_distance=parameters["facet_distance"],
    //                       CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
    //                       CGAL::parameters::cell_size=parameters["cell_size"] );



    //Mesh_criteria criteria(CGAL::parameters::cell_size=parameters["edge_size"],CGAL::parameters::facet_distance=parameters["facet_angle"]); // pointer ?? 

    Mesh_criteria criteria(CGAL::parameters::cell_size=parameters["cell_size"]);
    c3t3 = CGAL::make_mesh_3<C3t3>(*domain_ptr.get(), criteria); 

    //ERROR typedef CGAL::Mesh_3::internal::Edge_criteria_sizing_field_wrapper<Edge_criteria> Sizing_field;
    //CGAL::Mesh_3::Protect_edges_sizing_field<C3t3, Mesh_domain, Sizing_field>protect_edges(c3t3, *domain_ptr.get(), Sizing_field(criteria.edge_criteria_object()), 0.01);
    //protect_edges(true);

    remove_isolated_vertices(c3t3);

    


}
void CGALMeshCreator::create_mesh(int inital_points)
{
    Mesh_criteria criteria(CGAL::parameters::edge_size  =parameters["edge_size"],
                           CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );

    /* CGAL::internal::Mesh_3::init_c3t3(c3t3, *domain_ptr.get(), criteria, inital_points); */

    /* refine_mesh_3(c3t3, *domain_ptr.get(), criteria,CGAL::parameters::no_reset_c3t3()); */

    /* remove_isolated_vertices(c3t3); */

}

void CGALMeshCreator::save_mesh(std::string OutPath)
{
    std::ofstream  medit_file(OutPath);
    c3t3.output_to_medit(medit_file);
    medit_file.close();
}


void CGALMeshCreator::refine_mesh()
{

   Mesh_criteria criteria(CGAL::parameters::edge_size  =parameters["edge_size"],
                           CGAL::parameters::facet_angle=parameters["facet_angle"],
                           CGAL::parameters::facet_size =parameters["facet_size"],
                           CGAL::parameters::facet_distance=parameters["facet_distance"],
                           CGAL::parameters::cell_radius_edge_ratio=parameters["cell_radius_edge_ratio"],
                           CGAL::parameters::cell_size=parameters["cell_size"] );


   refine_mesh_3(c3t3, *domain_ptr.get(), criteria,CGAL::parameters::no_reset_c3t3());


}




void CGALMeshCreator::add_sharp_edges(Polyhedron& polyhedron) // rename store sharp edges -> cgal 4.13 try new command
{ 
  for( Polyhedron::Edge_iterator he = polyhedron.edges_begin(); he != polyhedron.edges_end() ; ++he)
  {

     if(he->is_feature_edge() )
     {      
        Polyline_3 polyline;
        polyline.push_back(he->vertex()->point());
        polyline.push_back(he->opposite()->vertex()->point());      
        add_polyline(polyline);
     }    
  }
}



template < typename InputIterator>
void CGALMeshCreator::insert_edge(InputIterator begin, InputIterator end) // Used in conjuction with refine mesh -> will insert points directly into mesh
{
  Tr& tr = c3t3.triangulation();
  Corner_index corner_index (1);
  Curve_index curve_index (1);
  std::vector<Vertex_handle> vertex_map;

  for ( ; begin != end; ++begin)
  {

      Vertex_handle *vh =  new Vertex_handle();
      *vh =  tr.insert(*begin);
      vertex_map.push_back(&vh);

  }
  std::cout << "Size: " <<  vertex_map.size() << '\n';

  for (typename std::vector<Vertex_handle>::const_iterator it = vertex_map.begin();it != vertex_map.end(); it++)
  {
      if (  std::next(it, 1) != vertex_map.end() )
      {
          c3t3.add_to_complex(*it,*std::next(it, 1), curve_index);  // split ??
          
      }

  } 


}

template < typename InputIterator> // Used in conjuection with make mesh, the 1d feature is presevered if possible in the meshing. -> Problem higer restriction causes errors 
void CGALMeshCreator::insert_edges(InputIterator begin, InputIterator end)
{

  for ( ; begin != end; ++begin)
  {
     CGALMeshCreator::add_polyline( begin->begin(), begin->end() );

  }
}


void CGALMeshCreator::label_boundary_cells(int btag , int ntag ) // workaround to mark boundary cells of for example lateral ventircles. this allows for easy marking of Facetfunction in FEniCS 
{
  Subdomain_index subdomain_index(btag);
  Subdomain_index subdomain_index_bis(ntag);
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin(subdomain_index);cit != c3t3.cells_in_complex_end(); ++cit)
  {
     for (std::size_t i = 0; i < 4; i++)
     { if (c3t3.subdomain_index(cit->neighbor(i))!=subdomain_index){c3t3.set_subdomain_index(cit->neighbor(i), subdomain_index_bis);}    }      
  }

} 

void CGALMeshCreator::remove_label_cells(int tag) //rename
{
  Subdomain_index subdomain_index(tag);
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin(subdomain_index);cit != c3t3.cells_in_complex_end(); ++cit){ c3t3.remove_from_complex(cit);}
}






