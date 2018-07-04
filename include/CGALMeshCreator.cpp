#include "CGALSurface.h"
#include "CGALMeshCreator.h"


int main() {
    auto foo = CGALSurface("hei");
    auto bar = CGALMeshCreator(foo);
    bar.default_parameters();

    std::cout << "Hello, world!" << std::endl;
}
