#ifndef  __Polyhedral_vector_to_labeled_function_wrapper_H

#define __Polyhedral_vector_to_labeled_function_wrapper_H



#include <boost/dynamic_bitset.hpp>


#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include "SubdomainMap.h"

#include "CGAL/Polyhedral_mesh_domain_3.h"

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

                //typedef class AABB_const_polyhedron_edge_primitive<BGT, Polyhedron_> AABB_primitive;
                //typedef class AABB_traits<BGT,AABB_primitive>     AABB_traits;
                //typedef class AABB_tree<AABB_traits>              AABB_tree;

                //tree_(new AABB_tree(p.edges_begin(), p.edges_end()))



                Polyhedral_vector_to_labeled_function_wrapper(std::vector<Function_*>& v) : function_vector_(v) 
                {
                    DefaultMap* map = new DefaultMap();
                    subdmap =map ;
                }
                
                Polyhedral_vector_to_labeled_function_wrapper(std::vector<Function_*>& v, AbstractMap& map) : function_vector_(v)
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

                //return_type index_from_surface_patch_index

                /*struct Construct_intersection
                {
                       Construct_intersection(const Polyhedral_vector_to_labeled_function_wrapper& domain) : r_domain_(domain) {std::cout<< "!he" << std::endl;}

                    private:
                        const Polyhedral_vector_to_labeled_function_wrapper& r_domain_;
                }; 

                Construct_intersection construct_intersection_object() const
                {
                    std::cout<< "!he" << std::endl;
                    return Construct_intersection(*this);
                }



                return_type operator()(const Segment_3& s) const
                { 
                   return this->operator()(s.source(), s.target());
                }
                return_type operator()(const Point_3& p1 , const Point_3& p2 ) 
                {
                    
                    // points do not match value 
                    int nb_func = function_vector_.size();
                    Bmask bits1(nb_func);
                    Bmask bits2(nb_func);
                    for ( int i = 0 ; i < nb_func ; ++i )
                    {
                        bits1[i] =(bool)function_vector_[i]->is_in_domain_object()(p1);
                        bits2[i] =(bool)function_vector_[i]->is_in_domain_object()(p2);
                    }
                    //std::cout << bits << std::endl;

                    std::cout<< " get" << std::endl;
                    //return subdmap->index(bits);
                    return static_cast<return_type>(bits1^bits2);
                }
                return_type surface_patch_index(const Index& index) 
                {
                    std::cout<< "!he" << std::endl;
                    return static_cast<return_type>(1);
                }*/

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
