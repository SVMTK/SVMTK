#ifndef __READ_POLYGONS_STL_H
#define __READ_POLYGONS_STL_H


#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include <math.h>


template<typename T>
double convert_string(const std::string& s)
{
    // TODO: inline
    std::istringstream is(s);
    T val;
    is >> val;

    return val;
}


void get_next_line(std::ifstream& file, std::string& line, std::size_t &lineno)
{
    // TODO: inline
    do
    {
        std::getline(file, line);
        boost::algorithm::trim(line);
        lineno++;
    }
    while (!file.eof() && line == "");
}


template<typename Point_3, typename  Polygon_3>
bool read_polygons_STL(std::ifstream& file,
            std::vector<Point_3>& points,
            std::vector<Polygon_3>& facets,
            bool verbose = false)
{
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    std::map<Point_3,int> pmap; 

    std::string facet("facet"),
        vertex("vertex"),
        endloop("endloop"),
        endsolid("endsolid");

    Polygon_3 tpolygon;
    Point_3 tpoint;

    std::string line;
    std::size_t lineno = 0;
    const boost::char_separator<char> sep(" ");

    get_next_line(file, line, lineno);
    if (line.substr(0, 5) != "solid")
    {
       std::cout<<" Error0" << line <<std::endl;
       return false;
    }

    get_next_line(file, line, lineno);

    int hva=0;
    int count = 1; // Shift 
    do
    {
        if ( line.substr(0, 5) == facet)
        {
            get_next_line(file, line, lineno);
            if (line != "outer loop")
            {
                std::cout<<" Error1" << line <<std::endl;
                return false;
            }

            get_next_line(file, line, lineno);
            tpolygon.clear();
            do
            {
                tokenizer tokens(line, sep);
                tokenizer::iterator tok_iter = tokens.begin();
                if (*tok_iter != "vertex")
                {
                    std::cout<<" Error2 " << line <<std::endl;
                    return false;
                }

                ++tok_iter;
                const double x = convert_string<double>(*tok_iter); ++tok_iter;
                const double y = convert_string<double>(*tok_iter); ++tok_iter;
                const double z = convert_string<double>(*tok_iter); ++tok_iter;

                tpoint = Point_3(x, y, z);
                if (pmap[tpoint] == 0) // if point is zero -> assign count
                {
                    pmap[tpoint] = count;
                    points.push_back(tpoint);
                    count++; 

                }

                tpolygon.push_back(pmap[tpoint] - 1);  // Reshift
                get_next_line(file, line, lineno);
            }
            while (line.substr(0,7)!=endloop); 

            get_next_line(file, line, lineno);
        }   

        facets.push_back(tpolygon);
        get_next_line(file, line, lineno);
    }
    while (line.substr(0, 8) != endsolid); 
  
    std::cout << facets.size() << std::endl;
    std::cout << points.size() << std::endl;
    std::cout << pmap.size() << std::endl;
    pmap.clear();
    return true;
}


#endif
