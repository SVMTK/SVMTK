#ifndef  __Polyhedral_vector_to_labeled_function_wrapper_H
#define __Polyhedral_vector_to_labeled_function_wrapper_H

#include <boost/dynamic_bitset.hpp>

#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include "SubdomainMap.h"

#include "CGAL/Polyhedral_mesh_domain_3.h"


namespace CGAL {

template<class Function_, class BGT>
class Polyhedral_vector_to_labeled_function_wrapper
{
    public:
        // Types
        typedef int return_type;  
        typedef std::vector<Function_*>   Function_vector; // Polyhedron_mesh_domain_with features -> stores polyhedron

        typedef typename BGT::Point_3       Point_3;
        typedef typename BGT::Segment_3 Segment_3;
        typedef typename Function_::Index Index;
        typedef boost::dynamic_bitset<>   Bmask;
        typedef typename BGT::FT FT;

        Polyhedral_vector_to_labeled_function_wrapper(std::vector< Function_* >& v): function_vector_(v)
        {
             DefaultMap* map = new DefaultMap();
             subdmap =map ;
        }

        Polyhedral_vector_to_labeled_function_wrapper(std::vector<Function_*>& v, AbstractMap& map) : function_vector_(v)
        {
            subdmap = j&map;
        }

        ~Polyhedral_vector_to_labeled_function_wrapper() {}

        return_type operator()(const Point_3& p, bool use_cache = false) const
        {
            int nb_func = function_vector_.size();
            Bmask bits(nb_func);

            for ( int i = 0 ; i < nb_func ; ++i )
            {
                bits[i] = (bool)function_vector_[i]->is_in_domain_object()(p);
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


}   // end namespace CGAL

#endif
