#ifndef Neuron_H


#define Neuron_H

#include "CGALSurface.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
struct NData // move into class ? 
{
    double x,y,z,radius;
    int index, next,type; 

    NData(){}

    ~NData(){}
    
    void print()
    {
       std::cout<< index << "\t" << type << "\t " << x << "\t" <<  y << "\t " << z << "\t" << radius<<"\t" <<next << std::endl; 
    }

    friend std::istream& operator>>(std::istream &is, NData&d)
    {
        return is >> d.index >> d.type >> d.x>> d.y>> d.z >> d.radius >> d.next;
    }

};


class Neuron : public CGALSurface
{

   // CGALSurface but with polylines and radius as input 

   // i.e MeshCreator(const Neuron& neuron) -> make mesh of composite cones -> require that MeshCreator gets_polyhedron/get_mesh inheritance


   // MeshCreator.add_polylines(const Neuron& neuron) -> add polylines -> possible to template, but get_polylines is also an option 

   // TODO: Add features to CGALSurface get set and clear mesh

   // Question :  should each segment be a CGALSurface ?

   public:

       typedef std::vector<NData>::iterator NData_iteator;
       typedef CGALSurface::Polylines Polylines;
       typedef CGALSurface::Polyline_3 Polyline_3;
       typedef CGALSurface::Point_3 Point_3;


       Neuron();

       ~Neuron() {}

       Neuron(std::string filename);

       Polylines get_features();
       
       void surface_segment(double xb,double yb,double zb, double xt,double yt,double zt,double rb,double rt); // struct input ?
         
       void surface_mesh();

       void print(){ for(std::vector<NData>::iterator sit= segments.begin();  sit!= segments.end(); ++sit) { sit->print();}}

       void num_lines(){ std::cout<< segments.size() << std::endl;}

   private :
       std::vector<NData> segments;
       
             


};
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
  //TDOImplemnt 

}
#endif
