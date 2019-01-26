#include "Neuron.h"

#include <CGAL/Polygon_mesh_processing/corefinement.h>
Neuron::Neuron(std::string filename)
{
    std::ifstream in(filename); //don't forget to check `if (in.is_open())` -> FIX
    std::vector<NData> temp{std::istream_iterator<NData>(in), std::istream_iterator<NData>()}; 

    segments = temp;
}

void Neuron::surface_mesh()
{

  std::vector<NData>::iterator nit; 
  int count=0;
  for(std::vector<NData>::iterator sit= std::next(segments.begin());  sit!= segments.end(); ++sit) // typedef iterator tree-structure first element not important | Maybe Sphere ?
  {
    
     // ADD if for type 
     nit =  std::next(segments.begin(),sit->next-1); 

     std::cout<< "####################################" << std::endl;
     sit->print();
     nit->print();
     std::cout<< "####################################" << std::endl;
     surface_segment(sit->x,sit->y,sit->z,nit->x, nit->y,nit->z, sit->radius, nit->radius);  



  }
} 

Neuron::Polylines Neuron::get_features()
{

  Polylines polylines;

  std::vector<NData>::iterator nit; 
  for(std::vector<NData>::iterator sit= std::next(segments.begin());  sit!= segments.end(); ++sit) 
  {
    Polyline_3 polyline;

    nit =  segments.begin() + sit->next; 
    
    polyline.push_back(Point_3(nit->x,nit->y,nit->z) );
    polyline.push_back(Point_3(sit->x,sit->y,sit->z) );
    polylines.push_back(polyline);
    
  }

  return polylines;
} 

void Neuron::surface_segment(double xb,double yb,double zb, double xt,double yt,double zt,double rb,double rt )// -> cast to CGALSURFACE
{
 
  CGALSurface body = CGALSurface();
  CGALSurface top = CGALSurface();

  body.make_cone(xb,yb,zb,xt,yt,zt,rb,rt,180);  // split.edges include for better points
  top.make_sphere(xb,yb,zb,rb);


  if (mesh.is_empty())
  {
    mesh = body.get_mesh(); 
   // mesh+=top.get_mesh();
  }
  else
  {
    mesh+=body.get_mesh(); 
    mesh+=top.get_mesh();
    //CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh, body.get_mesh(),mesh);
   // CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh, top.get_mesh(),mesh);
  }

}

