#ifndef  __Polyhedral_vector_to_labeled_function_wrapper_H
#define __Polyhedral_vector_to_labeled_function_wrapper_H


#include <boost/dynamic_bitset.hpp>


#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include "SubdomainMap.h"


namespace CGAL {
    template<class Function_, class BGT>
    class Polyhedral_vector_to_labeled_function_wrapper {
        public:
           // Types
           typedef int return_type;
           typedef std::vector<Function_*>   Function_vector;
           typedef typename BGT::Point_3     Point_3;
           /* typedef std::vector< boost::optional<int> > Bmask; */
           /* typedef std::vector<typename Function_::Subdomain_index>   Bmask; */
           /* typedef std::vector<int>   Bmask; */
           typedef boost::dynamic_bitset<>   Bmask;

            Polyhedral_vector_to_labeled_function_wrapper(std::vector<Function_*>& v) : function_vector_(v) {
                DefaultMap* map = new DefaultMap();
                subdmap =map ;
            }

            Polyhedral_vector_to_labeled_function_wrapper(
                    std::vector<Function_*>& v, const AbstractMap& map) : function_vector_(v) {
                subdmap = &map;
            }

            ~Polyhedral_vector_to_labeled_function_wrapper() {}

            return_type operator()(const Point_3& p, bool use_cache = false) const {
                int nb_func = function_vector_.size();
                Bmask bits(nb_func);
                boost::optional<int> tmp;
                for (int i = 0; i < nb_func; ++i) {
                    tmp = function_vector_[i]->is_in_domain_object()(p);
                    if (tmp) {
                       bits[i] = tmp.get();
                    }
                }

                return subdmap->index(bits);
            }

            Bbox_3 bbox() const {
                int nb_func = function_vector_.size();

                Bbox_3 sum_bbox;
                for (int i = 0 ; i < nb_func ; ++i) {
                    sum_bbox += function_vector_[i]->bbox();
                }
                return sum_bbox;
           }

        private:
            Function_vector function_vector_;
            AbstractMap* subdmap;
    };
}
#endif
