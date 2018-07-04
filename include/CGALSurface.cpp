#include "CGALSurface.h"


int main() {
    auto foo = CGALSurface("hei");
    foo.fill_holes();
    foo.triangulate_faces();
    foo.stitch_borders();
    foo.self_intersections();
    std::cout << "Hello, world!" << std::endl;
}
