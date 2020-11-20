#ifndef  __Polyhedral_vector_to_labeled_function_wrapper_H

#define __Polyhedral_vector_to_labeled_function_wrapper_H


#include "SubdomainMap.h" 
#include <CGAL/Polygon_mesh_processing/bbox.h>


/*
 * Wrapper : 
 * CGAL Polyhedron to CGAL labeled mesh.
 * Combines the option of using CGAL polyhedrons to create mesh 
 * with specific tags for overlapping surfaces. (CGAL version 4.11) 
 */
namespace CGAL {
        template<class Function_, class BGT>
        class Polyhedral_vector_to_labeled_function_wrapper
        {
            public:

                typedef int return_type;  
                typedef std::vector<Function_*>   Function_vector;
                typedef typename BGT::Point_3       Point_3;
                typedef boost::dynamic_bitset<>   Bmask;
            
                /**
                 * Constructor: 
                 * stores argument in member variables 
                 * @param v a vector of functions i.e. surfaces with query is inside  
                 * @param map a smart pointer to a child class of SVMTK virtuell class AbstractMap
                 *       options : DefaultMap and SubdomainMap
                 */
                Polyhedral_vector_to_labeled_function_wrapper(const std::vector<Function_*>& v, std::shared_ptr<AbstractMap> map) : function_vector_(v)
                {
                    subdmap =std::move(map);

                }

                ~Polyhedral_vector_to_labeled_function_wrapper() {}
                
               /**
                 * This function determins the tag of a point and subsequent enclosed cell
                 * by evaluating a bitstring linked to the vector of surfaces. The bitstring has a length equal to the number 
                 * of surfaces, and for each bit the value "0" means the point is outside 
                 * and "1" means inside the corresponding surface to bit position.
                 * This will ensure unqiueness for all combination of overlapping surfaces
                 * Then the tag is determined by a mapping og the bitstring based on the constructor argument map.
                 * @param p a CGAL 3D Point object 
                 * @param use_cache not in use 
                 * @return a mapping of the bistring to an integer based constructor argument map.
                 */
                return_type operator()(const Point_3& p, bool use_cache = false) const
                {
                    int nb_func = function_vector_.size();
                    Bmask bits(nb_func);

                    for ( int i = 0 ; i < nb_func ; ++i )
                    {
                        bits[i] =(bool)function_vector_[i]->is_in_domain_object()(p);
                    }

                    return subdmap->index(bits);
                }
                
               /**
                 * Constructs a bounding box of the surfaces.
                 * @return sum_bbox the summation of CGAL bounding boxes for each surfaces 
                 */
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
                std::shared_ptr<AbstractMap> subdmap;
        };
}



#endif
