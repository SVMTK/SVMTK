from svmtk import Surface, Domain


def chair_function(x, y, z):
    x2 = x*x
    y2 = y*y
    z2 = z*z
    x4 = x2*x2
    y4 = y2*y2
    z4 = z2*z2

    surf = x4 - 1.2*x2*y2 + 3.6*x2*z2 - 7.50*x2 + y4 + 3.6*y2*z2 - 7.50*y2 + .2*z4
    surf -= 7.50*z2 + 64.0625 - 16.0*z*y2 + 16.0*x2*z
    return surf


if __name__ == "__main__":
    surf = Surface()
    surf.implicit_surface(chair_function, 6.0, 30, 0.1, 0.1)
    surf.save("results/chair.off")

    maker = Domain(surf)
    maker.create_mesh()
    maker.save("results/chair.mesh")
