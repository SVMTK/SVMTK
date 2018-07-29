#include "CGALSurface.h"
#include "reconstruct_surface.h"


int main() {
    auto foo = CGALSurface("hei");
    foo.fill_holes();
    foo.triangulate_faces();
    foo.stitch_borders();
    foo.self_intersections();

    auto reconstruct_surface(foo);

    std::cout << "Hello, world!" << std::endl;
}
