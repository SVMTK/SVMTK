#include <catch.hpp>

#include "Domain.h"
#include "Surface.h"


TEST_CASE("Minimum bounding sphere")
{
    Surface surface; 
    surface.make_sphere(0.,0.,0.,3,20); 
    Domain domain(surface);
    REQUIRE( domain.get_bounding_sphere_radius()==Approx(3).margin(1e-3)) ; // less than 4 larger than 2 ?
}

