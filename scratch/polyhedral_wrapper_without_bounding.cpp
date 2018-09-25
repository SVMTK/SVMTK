#include <CGAL/license/Mesh_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>


#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include "make_multicomponent_mesh_3.h"

#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Mesh_3/Dump_c3t3.h>
#include <CGAL/IO/File_medit.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>

#include <functional>

#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/config.h>
#include <CGAL/assertions.h>
// Domain Mesh_polyhedron_3
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Polyhedral_mesh_domain_3;
typedef boost::dynamic_bitset<>   Bmask;

class AbstractMap
{
   public:
        typedef int return_type;
        typedef boost::dynamic_bitset<> Bmask;
        virtual return_type index(const Bmask bits) = 0;
        
        AbstractMap() {};
        ~AbstractMap() {}; 
};


class DefaultMap : public AbstractMap
{
   public:
        typedef int return_type;
        typedef boost::dynamic_bitset<> Bmask;
    
        DefaultMap() {}
        ~DefaultMap() {} 

        return_type index(const Bmask bits) 
        {
           return static_cast<return_type>(bits.to_ulong());
        }
};

class SubdomainMap : public AbstractMap
{
   private:
        std::map<boost::dynamic_bitset<>,int> subdmap;

   public:
        typedef int return_type;
        typedef boost::dynamic_bitset<> Bmask;
    
        SubdomainMap() {}
        ~SubdomainMap() {} 

        void add(std::string string, int subdomain)
        {
           subdmap[Bmask(string)]=subdomain;
        } 
        return_type index(const Bmask bits) 
        {
           if (subdmap.find(bits)!=subdmap.end()){
              return static_cast<return_type>(subdmap[bits]); 
           } 
           else{
              
              return static_cast<return_type>(0); 
           }

           
        }
        void print() 
        {
           for(auto elem : subdmap)
           {
              std::cout << elem.first << " " << elem.second << " " << std::endl;
           }
        }
};



namespace CGAL {
        template<class Function_, class BGT>
        class Polyhedral_vector_to_labeled_function_wrapper
        {
            public:
                // Types
                typedef int return_type;
                typedef std::vector<Function_*>   Function_vector;
                typedef typename BGT::Point_3	  Point_3;
                typedef boost::dynamic_bitset<>   Bmask;


                typedef std::function<return_type(Bmask)> fbind;
                typedef  std::map<boost::dynamic_bitset<>,int> domainmap;



                Polyhedral_vector_to_labeled_function_wrapper(std::vector<Function_*>& v) : function_vector_(v)  // , generator()
                {
                    subdmap = new DefaultMap();

                }

                Polyhedral_vector_to_labeled_function_wrapper(std::vector<Function_*>& v, AbstractMap& map) : function_vector_(v) //fbind (fptr*) (fbind*,return_type)
                {
                    subdmap =&map;

                }
                ~Polyhedral_vector_to_labeled_function_wrapper() {}

                return_type operator()(const Point_3& p, bool use_cache = false) const
                {
                    int nb_func = function_vector_.size();
                    Bmask bits(nb_func);

                    for ( int i = 0 ; i < nb_func ; ++i )
                    {
                        bits[i]=(bool)function_vector_[i]->is_in_domain_object()(p); // change ++
                    }

                    //std::cout << subdmap->index(bits) << bits << std::endl;
                    return subdmap->index(bits);
                }

                Bbox_3 bbox() const
                {
                    int nb_func = function_vector_.size();

                    Bbox_3 sum_bbox;
                    for ( int i = 0 ; i < nb_func ; ++i )
                    {
                        sum_bbox+= function_vector_[i]->bbox();
                    }
                    return sum_bbox;
                }


            private:
                Function_vector function_vector_;
                fbind domain_index;
                AbstractMap* subdmap;
        };
}




#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, K> Function_wrapper;// NAMESPACE
typedef Function_wrapper::Function_vector Function_vector; 
typedef CGAL::Labeled_mesh_domain_3<Function_wrapper, K> Mesh_domain; 

// Triangulation

typedef CGAL::Mesh_triangulation_3<Polyhedral_mesh_domain_3,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Tr::Vertex_handle Vertex_handle;
typedef Tr::Cell_handle Cell_handle;
using namespace CGAL::parameters;

typedef C3t3::Triangulation::Facet Facet; 
typedef C3t3::Subdomain_index Subdomain_index;
typedef C3t3::Surface_patch_index Surface_patch_index;
int main(int argc, char*argv[])
{
  const char* fname;
  const char* OutPath;
  const char* fname_b=argv[1];

  Function_vector v ; 

  for (int i = 1; i < argc-1; i++) 
  { 
     fname = argv[i];
     std::ifstream input(fname);
     Polyhedron *polyhedron = new Polyhedron();
     input >> *polyhedron;
     if(input.fail())
     {
       std::cerr << "Error: Cannot read file" <<  fname << std::endl;
       return EXIT_FAILURE;
     }
     input.close();
     std::cout << polyhedron << std::endl;

     Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(*polyhedron);

     v.push_back(polyhedral_domain);
  
  }

  

  OutPath = argv[argc-1];

  std::string sstring1("0001");

  std::string sstring2("0011");
  std::string sstring3("0101");
  std::string sstring4("0111");
  std::string sstring5("0100");

  std::string sstring6("1011");
  std::string sstring7("1101");
  std::string sstring8("1001");
  std::string sstring9("1000");
  //std::string sstring5("1100");
  
  //std::string sstring4("000111");
  //std::string sstring6("111000");
  //std::map<boost::dynamic_bitset<>,int> ;

  SubdomainMap subdomain_map;


  subdomain_map.add(sstring2,2);
  subdomain_map.add(sstring3,3);
  subdomain_map.add(sstring4,3);

  subdomain_map.print();

  std::vector<std::string> vps; 


  //vps.push_back(sstring4);
  //vps.push_back(sstring5);
  //vps.push_back(sstring6);
  

  Function_wrapper wrapper = Function_wrapper(v,subdomain_map);


  Mesh_domain domain(wrapper,wrapper.bbox() ); // ?? 
  double cellsize= 1.8;
                                       
  Mesh_criteria criteria(edge_size = cellsize , facet_angle=25, facet_size=cellsize, facet_distance=cellsize/10,cell_radius_edge_ratio=3,cell_size=cellsize );

  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria ); 


  // TODO: Error making mesh where number of points is not sufficient -> polyhedron should have as many points as possible -> subdiv -> until 1000

  //C3t3 c3t3;
  //make_multicomponent_mesh_3_impl<C3t3>(c3t3,domain,criteria,CGAL::parameters::no_exude(),CGAL::parameters::no_perturb(),CGAL::parameters::no_odt(),CGAL::parameters::no_lloyd(),false);

  Subdomain_index subdomain_index(7);
  Subdomain_index subdomain_index_bis(8);
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin(subdomain_index);cit != c3t3.cells_in_complex_end(); ++cit)
  {

     for (std::size_t i = 0; i < 4; i++)
      {   
          if (c3t3.subdomain_index(cit->neighbor(i))!=subdomain_index){
             c3t3.set_subdomain_index(cit->neighbor(i), subdomain_index_bis);
          }    
      }
      
  }


  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin(subdomain_index);cit != c3t3.cells_in_complex_end(); ++cit)
  {

      c3t3.remove_from_complex(cit);
  }

  std::ofstream outputFacets("Facet.txt",std::ofstream::out);


  int counter=0;
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
      counter= counter+ 1;
      for (std::size_t i = 0; i < 4; i++)
      {   
          if (cit->surface_patch_index(i).first==7 or cit->surface_patch_index(i).second==7){
          outputFacets <<  counter << " "<<  i <<" "<< 14 << std::endl;
          }
          else {
          outputFacets <<  counter << " "<<  i <<" "<< 0 << std::endl;
          }
      }
  }
  outputFacets.close();



  // The following code can be optimized
  std::map<Vertex_handle, bool> vertex_map;
  for( C3t3::Triangulation::Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin();vit != c3t3.triangulation().finite_vertices_end();++vit)
  { 
      vertex_map[vit] = false ;  
  }
 
  for(C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
    for (std::size_t i = 0; i < 4; i++)
    {
        vertex_map[cit->vertex(i)] = true;
    }
  }

  std::cout<<  "Num vertices" << c3t3.triangulation().number_of_vertices ()  << std::endl;
  
  for (std::map<Vertex_handle, bool>::const_iterator it = vertex_map.begin();it != vertex_map.end(); it++)
  {
    if (!it->second) 
    {
       c3t3.triangulation().remove(it->first);
    }
  }

  
  std::cout<< "Num vvertices " <<  c3t3.triangulation().number_of_vertices ()  << std::endl;

  std::ofstream off_file("out.off");
  c3t3.output_boundary_to_off(off_file);
  std::ofstream *medit_file= new std::ofstream(OutPath);
  std::cout <<c3t3.number_of_cells() << std::endl; 
  output_to_medit(*medit_file,c3t3,true,true);
  medit_file->close();

  return 0;


}
