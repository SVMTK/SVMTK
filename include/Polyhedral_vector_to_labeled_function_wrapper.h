#ifndef  __Polyhedral_vector_to_labeled_function_wrapper_H

#define __Polyhedral_vector_to_labeled_function_wrapper_H


#include "SubdomainMap.h" // is comipled before this safegaurd capture 
//s#include <boost/dynamic_bitset.hpp>


//#include <list>
//#include <string>
//#include <iostream>
//#include <fstream>


//#include "CGAL/Polyhedral_mesh_domain_3.h" ??

//#include <CGAL/AABB_tree.h>
//#include <CGAL/AABB_traits.h>
//#include <CGAL/AABB_polyhedron_triangle_primitive.h>
//#include <CGAL/AABB_polyhedron_segment_primitive.h>

//#include <CGAL/point_generators_3.h>



namespace CGAL {
        template<class Function_, class BGT>
        class Polyhedral_vector_to_labeled_function_wrapper
        {
            public:
                // Types
                typedef int return_type;  

                typedef std::vector<Function_*>   Function_vector; // Polyhedron_mesh_domain_with features -> stores polyhedron
                // TODO
                // Polyhedron_mesh_domain_with features  operator() 
                // similar to Hybrid mesh domain
                // SubID generator 
         
                typedef typename BGT::Point_3       Point_3;
                typedef typename BGT::Segment_3 Segment_3;
                typedef typename Function_::Index Index;
                typedef boost::dynamic_bitset<>   Bmask;
                //typedef typename Function_::FT    FT;
                typedef typename BGT::FT FT;
                 

                Polyhedral_vector_to_labeled_function_wrapper(const std::vector<Function_*>& v, AbstractMap& map) : function_vector_(v)
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
                        bits[i] =(bool)function_vector_[i]->is_in_domain_object()(p);
                    }

                    return subdmap->index(bits);
                }
                return_type operator()(const Segment_3& segment) const
                {
                    int nb_func = function_vector_.size();
                    Bmask bits(nb_func);
                    std::cout << " XXX " << std::endl; 
                    for ( int i = 0 ; i < nb_func ; ++i )
                    {
                        bits[i] =(bool)function_vector_[i]->do_intersect_surface_object()(segment);
                    }

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
                AbstractMap* subdmap;
        };
}



#endif
