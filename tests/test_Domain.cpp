#include <catch.hpp>

#include "Surface.h"
#include "Domain.h"



TEST_CASE("Minimum bounding sphere")
{
    Surface surface; 
    surface.make_sphere<Domain>(0.,0.,0.,3,0.5,1.0e-6); 
    Domain domain(surface);
    REQUIRE( domain.get_bounding_sphere_radius()==Approx(3).margin(1e-3)) ; // less than 4 larger than 2 ?
}

