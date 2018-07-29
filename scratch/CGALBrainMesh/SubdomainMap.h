#ifndef SubdomainMap_H


#define SubdomainMap_H

#include <boost/dynamic_bitset.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_cv.hpp>
#include <map>
#include <string>
#include <vector>
#include <iostream>

class AbstractMap
{
   public:
        typedef int return_type;
        typedef boost::dynamic_bitset<> Bmask;
        virtual return_type index(const Bmask bits) = 0;
        
        AbstractMap() {}
        ~AbstractMap() {}
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
           return static_cast<return_type>(subdmap[bits]);  
        }
        void print() 
        {
           for(std::map<boost::dynamic_bitset<>,int>::iterator it=subdmap.begin();it!=subdmap.end();++it )
           {
              std::cout << it->first << " " << it->second << " " << std::endl;
           }
        }
};

#endif
