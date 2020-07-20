// Copyright (C) 2018-2019 Lars Magnus Valnes and 
//
// This file is part of Surface Volume Meshing Toolkit (SVM-TK).
//
// SVM-Tk is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SVM-Tk is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SVM-Tk.  If not, see <http://www.gnu.org/licenses/>.
#ifndef SubdomainMap_H


#define SubdomainMap_H

#include <boost/dynamic_bitset.hpp>
//#include <boost/type_traits/remove_reference.hpp> // ?
//#include <boost/type_traits/remove_cv.hpp> // ? 
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
        virtual ~AbstractMap() {}

};


class DefaultMap : virtual public AbstractMap
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

class SubdomainMap :virtual public AbstractMap
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
           std::reverse( string.begin(), string.end() ) ; // Bmask initiate reverse 
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
        int number_of_domains()
        {
           return subdmap.size() ; 


        }
      
};



#endif
