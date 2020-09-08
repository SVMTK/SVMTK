#ifndef  __Polyhedral_vector_to_labeled_function_wrapper_H

#define __Polyhedral_vector_to_labeled_function_wrapper_H


#include "SubdomainMap.h" 
#include <CGAL/Polygon_mesh_processing/bbox.h>

namespace CGAL {
        template<class Function_, class BGT>
        class Polyhedral_vector_to_labeled_function_wrapper
        {
            public:

                typedef int return_type;  
                typedef std::vector<Function_*>   Function_vector;
                typedef typename BGT::Point_3       Point_3;
                typedef boost::dynamic_bitset<>   Bmask;
            
                Polyhedral_vector_to_labeled_function_wrapper(const std::vector<Function_*>& v, std::shared_ptr<AbstractMap> map) : function_vector_(v)
                {
                    subdmap =std::move(map);

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
