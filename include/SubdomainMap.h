// Copyright (C) 2018-2021 Lars Magnus Valnes 
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
/* --- Includes -- */
#include "Errors.h"

/* -- STL -- */
#include <map>
#include <string>
#include <vector>
#include <iostream>

/* -- boost-- */
#include <boost/dynamic_bitset.hpp>
#include <boost/lexical_cast.hpp>




/**
 * \class 
 * The Abstract superclass for SubdomainMap.
 */
class AbstractMap
{
   public:
        typedef int return_type;
        typedef boost::dynamic_bitset<> Bmask;
        
        virtual return_type index(const Bmask bits) = 0;
        virtual const std::map<std::pair<int,int>,int> make_interfaces(std::vector<std::pair<int,int>> interfaces)=0;
        AbstractMap() {}
        virtual ~AbstractMap() {}

};

/**
 *  \class 
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
         * @brief Maps bistring to an integer using binary conversion to integer
         * @param bits a bitstring of the form boost::dynamic_bitset.
         * @return long conversion of the bit string.
         */
        return_type index(const Bmask bits) 
        {
           return static_cast<return_type>(bits.to_ulong());
        }
        struct sort_pairs{
                  template<typename T>
                  bool operator()(const std::vector<T> & a, const std::vector<T> & b)
                       { return a.size() > b.size(); }                   
       };           
        /** 
         * @brief Gives each interfaces an unique tag based on presence in mesh and 
         * from the highest cell tag. Ex:(Cell tags 1 2 -> Interface tags 3 4 5)
         * @param interfaces a vector of integer pairs that represent cell tag between a facet.
         * @return patches a map of unique combination of pairs with an unique integer. 
         */      
        const std::map<std::pair<int,int>,int> make_interfaces(std::vector<std::pair<int,int>> interfaces)
        { 
        
           std::map<std::pair<int,int> ,int> patches;  
           int iter =0;
           for( auto interface : interfaces )  
           {  
               if (interface.first > iter)
                    iter = interface.first ;                 
           } 
           iter++;
           for( auto interface : interfaces )  
           {  
               patches[std::pair<int,int>(interface.first,interface.second)]=iter++;
           } 
           return patches;
        }
};

/**
 * @brief Creates all possible binarystring given the number of elements elements
 * @param[in] string unfinished binarystirng, i.e. elements < max_elements 
 * @param[out] strings contains all possible binarystirng given the number of elements in binary string.
 * @param[in] max_elements number of elements in bitstring
 * @return void.
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
/**
 * \class
 *
 *  User specific method to set subdomains in the mesh.
 */
class SubdomainMap :virtual public AbstractMap
{

   public:
        typedef int return_type;
        typedef boost::dynamic_bitset<> Bmask;
       
        SubdomainMap(int number_of_surfaces) : num_surfaces(number_of_surfaces) 
        {
            //if (number_of_surfaces<1)
            //   throw InvalidArgumentError("Invalid number of surfaces"); 
            std::string zero_tag = std::string(number_of_surfaces , '0');
            add(zero_tag,0);        
        }
        //SubdomainMap() : num_surfaces(0) {}
         ~SubdomainMap() {} 

        /**
         * @brief Sets the number of surfaces
         *
         * The number of surfaces is used to automatically fill binary 
         * combinations when adding a sting with asterix.
         * @note Removing this function may cause errors in backwards compatiblity.
         * 
         * @param number_of_surfaces integer that equal the number of surfaces.
         * @return none.
         */
        void set_number_of_surfaces(int number_of_surfaces) 
        {
            if (number_of_surfaces<1)
               throw InvalidArgumentError("Invalid number of surfaces"); 
            this->num_surfaces= number_of_surfaces;
            std::string zero_tag = std::string(number_of_surfaces , '0');
            add(zero_tag,0);
        }

        /**
         * @brief Automatically fills the binary combintation of a string
         * with asterix *.
         * 
         * Precondition is that the number of surfaces is set to be non-zero.
         * Fils the remaining binary possibilites given a substring with asterix:  
         * i.e. 2 surfaces and *1 -> 01 and 11 
         *      3 surfaces and 11* -> 110 and 111 
         * The binary string containing only zero is set to zero in all cases. 
         * 
         * @param istring substring without * 
         * @param pos the position of * in substring before removal. 
         * @param tag the subdomain tag that is used for the substring 
         * @return none.
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
         * TODO FIXME
         * @brief Transforms a string of 0 and 1 to a binarystring 
         * 
         * @param string that contains only 0 and 1. 
         * @param subdomain an integer that determines the output tag of the subdomain
         * @return none 
         */
        void add(std::string string, int subdomain)
        {
        
           if ( string.find("*") !=std::string::npos)
           {
              if (this->num_surfaces==0)
                 throw InvalidArgumentError("Use of asterix require that number of surfaces to set.");    
           
              int pos = static_cast<int>(string.find("*"));
              string.erase(pos,1);
               // TODO better code writing 
              /*if (pos==0)

               else
                  string.erase(pos);
                */    
              if( static_cast<int>( string.length())<this->num_surfaces)
              {                                            
                  fill(string,pos,subdomain); 
              }
              else 
                 throw InvalidArgumentError("Bitstring exceeds number of surfaces");                   
           }
           else 
           {
              if(static_cast<int>( string.length())==this->num_surfaces or this->num_surfaces==0) 
              { 
                 std::reverse( string.begin(), string.end());
                 if (subdmap.find(Bmask(string)) == subdmap.end() )
                     subdmap[Bmask(string)]=subdomain;   
              }
              else 
                throw InvalidArgumentError( "Use the correct number of characters in bitstring." ); 
           }
        }      
        /**
         * @brief Erase binary string from SubdomainMap 
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
         * @brief Checks wether a point is inside a series of surfaces and returns 
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
         * @brief Prints Subdomains and Patches
         * @param none 
         * @return void
         *
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
         * @brief Returns all tags that is added to the class object.  
         * @param none 
         * @return tags vector of integer that represents the tags added to the class object.
         */
        std::vector<int> get_tags() 
        {
              std::vector<int> tags;
              for (auto it : subdmap) 
              {
                 tags.push_back(it.second);
              } 
              return tags;
        }
        
        /** 
         * @brief Returns the map. 
         * @param none 
         * @return tags vector of integer that represents the tags added to the class object.
         */
        std::map<std::string,int> get_map() 
        {
           std::map<std::string,int> result; 
           for(std::map<boost::dynamic_bitset<>,int>::iterator it=subdmap.begin();it!=subdmap.end();++it )
           {
              // Reversed so to look like what was added
              std::string dummy = boost::lexical_cast<std::string>(it->first);
              std::reverse(dummy.begin(),dummy.end());
              result[dummy] = it->second;
           }  
           return result;  
        }

        /**
         * @brief Adds a tag value for surfaces patches between subdomains defined by a pair of integer 
         * the class member variable patches 
         *
         * @param interface a integer pair that represents the suface interface between two subdomains 
         * @param tag the value used to represent the surface interface in the save file. 
         * @return void
         */
        void add_interface(std::pair<int,int> interface, int tag)
        {
              if (interface.second > interface.first) 
              {std::swap(interface.first, interface.second);}
              patches[interface] = tag;
        } 

        /**    
         * @brief Returns the content of class member variable patches between subdomains with 
         * the corresponding tag value. If patches is empty, gives each interfaces an unique 
         * tag based on presence in mesh and the highest cell tag 
         * @param interfaces contains a set of all interfaces between cells used in the mesh. 
         * @return patches a map  with pair of integers as key and a integer tag value.
         */
        const std::map<std::pair<int,int>, int> make_interfaces(std::vector<std::pair<int,int>> interfaces)
        {
           if (!patches.empty()) 
           {
              return patches;
           }
           else
           {  
              //FIXME
              
              int iter =0;
              for( auto interface : interfaces )  
              {  

                  if (interface.first > iter)
                      iter = interface.first;  
              } 
              /*auto tags =get_tags();
              auto iter = *std::max_element(std::begin(tags), std::end(tags));   */    
              iter++;

            
              for( auto interface : interfaces )  
              {  
                  if (patches.find(std::pair<int,int>(interface.first,interface.second)) == patches.end() )
                  {
                      patches[std::pair<int,int>(interface.first,interface.second)]=iter++;       
                  }    
              } 
              return patches;
           } 
        }
        
        /**
         * @brief get stored interface patches and tags.
         * @param none 
         * @param  
         */ 
        std::map<std::pair<int,int>, int> get_interfaces()
        {
           return patches;
        }

   private:
        int num_surfaces; 
        std::map<boost::dynamic_bitset<>,int> subdmap;
   protected:
        std::map<std::pair<int,int> ,int> patches;
};

#endif
