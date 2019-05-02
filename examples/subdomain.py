from svmtk import Surface, Domain


if __name__ == "__main__":
    sphere = Surface()
    sphere.make_cylinder(0., 0., 0., 2., 4., 3., 1.0, 21)

    sphere.split_edges(0.3)
    maker = Domain(sphere)
    maker.add_sharp_border_edges(sphere)
    maker.set_borders()
    maker.default_creating_mesh()
    maker.save("results/cylinder.mesh")
