  



#include <map>
#include <vector> 

// TODO: move
template<typename C3T3>
void remove_isolated_vertices(C3T3& c3t3)
{ 

  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;


  std::map<Vertex_handle, bool> vertex_map;
  for( Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin();vit != c3t3.triangulation().finite_vertices_end();++vit)
  { 
      vertex_map[vit] = false ;  
  }
 
  for(Cell_iterator cit = c3t3.cells_in_complex_begin();cit != c3t3.cells_in_complex_end(); ++cit)
  {
    for (std::size_t i = 0; i < 4; i++)
    {
        vertex_map[cit->vertex(i)] = true;
    }
  }

  int before = c3t3.triangulation().number_of_vertices() ;
  
  for (typename std::map<Vertex_handle, bool>::const_iterator it = vertex_map.begin();it != vertex_map.end(); ++it) // check post or pre increment
  {
    if (!it->second) 
    {
       c3t3.triangulation().remove(it->first);
    }
  }
  int after = c3t3.triangulation().number_of_vertices() ; 
  std::cout<<"Number of vertices removed: "  << before - after  << std::endl;

}
