from svmtk import Surface, Domain, Point_3


if __name__ == "__main__":
    surf = Surface()
    surf.make_cube(-2., -2., -2., 2., 2., 2.)

    maker = Domain(surf)
    maker.add_sharp_border_edges(surf)

    line1 = [Point_3(0, 0, -1.0), Point_3(0, 0, 0.0), Point_3(0, 1.0, 1.0)]
    line2 = [Point_3(0, 0, 0.0), Point_3(0, -1.0, 1.0)]

    maker.add_feature(line1)
    maker.add_feature(line2)

    maker.set_borders()
    maker.set_features()

    maker.create_mesh(24.0)
    maker.excude(100, 0)
    maker.save("bifurcation_in_cube.mesh")
