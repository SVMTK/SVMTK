// Copyright (C) 2018-2019 Lars Magnus Valnes and 
//
// This file is part of Surface Volume Meshing Toolkit (SVM-TK).
//
// SVM-Tk is free software: you can redistribute it and/or modify
// it undXFV4GHXgEzjYPm9er the terms of the GNU General Public License as published by
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
#include <boost/lexical_cast.hpp>
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

/**
 * Class: 
 * The Abstract superclass for SubdomainMap.
 */
class AbstractMap
{
   public:
        typedef int return_type;
        typedef boost::dynamic_bitset<> Bmask;
        
        virtual return_type index(const Bmask bits) = 0;
        //virtual return_type patch_index(const double s1,const double s2)=0;
        virtual const std::map<std::pair<int,int>,int> get_interfaces(const int number_of_surfaces)=0;
        //AbstractMap(int num_surfaces) :  _num_surfaces(num_surfaces)   {}
        AbstractMap() {}
        //get_num_surfaces() {return _num_surfaces;}
        virtual ~AbstractMap() {}

};

/**
 *  Class: 
 *  The defualt method to set subdomains in the mesh.
 *  Uses bitstring to integer conversion to set subdomain tag.
 *  
 */
class DefaultMap : virtual public AbstractMap
{
   public:
        typedef int return_type;
        typedef boost::dynamic_bitset<> Bmask;
    
        DefaultMap()  {}
        ~DefaultMap() {} 
         
        /** 
         * Maps bistring to an integer using binary conversion to integer
         * @param bits a bitstring of the form boost::dynamic_bitset.
         * @return long conversion of the bit string.
         */
        return_type index(const Bmask bits) 
        {
           return static_cast<return_type>(bits.to_ulong());
        }
              
        /** 
         * Return all possbile combination of a binarystring with number of elements  
         * equal to the number of surfaces.
         * @param integer of the number of surfaces 
         * @return patches a map of unique combination of pairs with an unique integer. 
         */      
        const std::map<std::pair<int,int>,int> get_interfaces(const int number_of_surfaces)
        { 
             
           int ulim =  pow (2,number_of_surfaces )+1; 
           std::map<std::pair<int,int> ,int> patches;
           int iter=1;
           for(int i =1; i< ulim; i++ )
           {
              for(int j=0; j< i ; j++)
              {
                 patches[std::pair<int,int>(i,j)]=iter++;
              }
           } 
           return patches;
        }
};

/**
 * Creates all possible binarystring given the number of elements elements
 * @param[in] string unfinished binarystirng, i.e. elements < max_elements 
 * @param[out] strings contains all possible binarystirng given the number of elements in binary string.
 * @param[in] max_elements number of elements in bitstring
 */
inline 
void generate_binary_strings(std::string string, std::vector<std::string>& strings, int max_elements) 
{
    if( static_cast<int>(string.length())<max_elements)
    {
       generate_binary_strings(string+"0", strings, max_elements); 
       generate_binary_strings(string+"1", strings, max_elements); 
    }
    else if ( static_cast<int>(string.length())==max_elements) 
    {
       strings.push_back(string);
    }    
    return;
}

class SubdomainMap :virtual public AbstractMap
{

   public:
        typedef int return_type;
        typedef boost::dynamic_bitset<> Bmask;
       
        SubdomainMap(int number_of_surfaces) : num_surfaces(number_of_surfaces) 
        {
            std::string zero_tag = std::string(number_of_surfaces , '0');
            add(zero_tag,0);        
        }
        SubdomainMap() : num_surfaces(0) {}
         ~SubdomainMap() {} 

        /**
         * Set the number of surfaces, optional 
         * @param number_of_surfaces integer that equal the number of surfaces.
         * @return none.
         */
        void set_number_of_surfaces(int number_of_surfaces) 
        {
            this->num_surfaces= number_of_surfaces;
            std::string zero_tag = std::string(number_of_surfaces , '0');
            add(zero_tag,0);
        }

        /**
         * @req that the number of surfaces is set to be non-zero.
         * Fils the remaining binary possibilites given a substring with the char *  
         * i.e. 2 surfaces and *1 -> 01 and 11 
         *      3 surfaces and 11* -> 110 and 111 
         * 
         * 
         * @param istring substring without * 
         * @param pos the position of * in substring before removal. 
         * @param tag the subdomain tag that is used for the substring 
         * @return none.
         * @note the binary string containing only zero is set to zero in all cases. 
         */
        void fill(std::string istring, int pos, int tag)
        {
             std::vector<std::string>  strings; 
             std::string dummy;
             generate_binary_strings(dummy, strings, this->num_surfaces - static_cast<int>(istring.length())  );  
             for ( auto i : strings)
             {   
                if (pos==0)
                { 
                   this->add(i+istring,tag); 
                }
                else 
                {  
                   this->add(istring+i,tag); 
                }
             }
        }

        /**
         * 
         * Transforms a string of 0 and 1 to a binqrystring 
         * 
         * @param string that contains only 0 and 1. 
         * @param subdomain an integer that determines the output tag of the subdomain
         * @return none 
         */
        void add(std::string string, int subdomain)
        {
           if( static_cast<int>( string.length()) !=this->num_surfaces and this->num_surfaces!=0 )
           {           
              if ( string.find("*") !=std::string::npos)
              {
                  int pos = static_cast<int>(string.find("*"));  
                  if (pos==0)
                     string.erase(pos,1);
                  else
                     string.erase(pos);                     
                  fill(string,pos,subdomain);
              }
              else 
              {
                 throw std::invalid_argument( "Number of surfaces do not match input string" );
              }              
           }  
           else 
           {
           // reverse the string since boost dyamic bitset reads the string in reverse
           std::reverse( string.begin(), string.end());
           if (subdmap.find(Bmask(string)) == subdmap.end() )
               subdmap[Bmask(string)]=subdomain;
           }
        } 
        
        /**
         * Erase binary string from SubdomainMap 
         * @param string a bitstring to remove from SubdomainMap 
         * @return void
         */
        void erase(std::string string)
        {
           // reverse the string since boost dyamic bitset reads the string in reverse
           std::reverse( string.begin(), string.end());
           subdmap.erase(Bmask(string));
        }
        /** 
         * Checks wether a point is inside a series of surfaces and returns 
         * a bitstring indicating with 1 if point is inside surface or 0 if not.
         * The bitstirng is used as a key in map that returns the added value for that 
         * specific region.
         * @param bits a bitstring of the form boost::dynamic_bitset.
         * @return return_type integer determined by a map with Bmask keys.
         */
        return_type index(const Bmask bits) 
        {
           return static_cast<return_type>(subdmap[bits]);  
        }
        /** 
         * Prints Subdomains and Patches
         */
        void print() 
        {
           for(std::map<boost::dynamic_bitset<>,int>::iterator it=subdmap.begin();it!=subdmap.end();++it )
           {
              // Reversed so to look like what was added
              std::string dummy = boost::lexical_cast<std::string>(it->first);
              std::reverse(dummy.begin(),dummy.end());
              std::cout<< "Subdomain: " << dummy << " " << it->second << " " << std::endl;
           }
           for(std::map<std::pair<int,int>,int>::iterator it=patches.begin();it!=patches.end();++it )
           {
              std::cout<< "Patches: " << it->first.first << " " << it->first.second << " " << it->second << " " << std::endl;
           }

        }
        /** 
         * Returns all tags that is added to the class object.  
         * @return tags vector of integer that represents the tags added to the class object.
         */
        std::vector<int> get_tags() 
        {
              std::vector<int> tags;
              tags.push_back(0);
              for (auto it : subdmap) 
              {
                 tags.push_back(it.second);
              } 
              return tags;
        }

        /**
         * Adds a tag value for surfaces patches between subdomains defined by a pair of integer 
         * the class member variable patches 
         * @param interface a integer pair that represents the suface interface between two subdomains 
         * @param tag the value used to represent the surface interface in the save file. 
         */
        void add_interface(std::pair<int,int> interface, int tag)
        {
              if (interface.second > interface.first) 
              {std::swap(interface.first, interface.second);}
              patches[interface] = tag;
        } 

        /**    
         * Returns the content of class member variable patches between subdomains with 
         * the corresponding tag value. If patches is empty, it will return all combinations 
         * of surfaces patches dependent on the number of surfaces.
         * @param number_of_surfaces an integer that indicates the number of surfaces 
         * @return patches a map  with pair of integers as key and a integer tag value.
         */
        const std::map<std::pair<int,int>, int> get_interfaces(const int number_of_surfaces)
        {
           if (!patches.empty()) 
           {
              return patches;
           }
           else
           {
              int iter=1;
              std::vector<int> tags = get_tags();
              for( auto i : tags )
              {
                 for(auto j : tags)
                 {
                    if( j>i )
                    {
                      patches[std::pair<int,int>(j,i)]=iter++;
                    }
                 } 
              }
              return patches;
           } 
        }

   private:
        int num_surfaces; 
        std::map<boost::dynamic_bitset<>,int> subdmap;
   protected:
        std::map<std::pair<int,int> ,int> patches;
};

#endif
