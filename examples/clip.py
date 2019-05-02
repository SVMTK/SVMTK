from svmtk import Surface, Domain


def pvs(x, y, z):
    return x**2 + 4.0*y**2 + x**2*y**2 - 4.0


if __name__ == "__main__":
    surf = Surface()
    surf.implicit_surface(pvs, 6.0, 30, 0.1, 0.1)
    surf.clip(0, 0, -1, -3, False)
    surf.clip(0, 0, 1, -3, False)
    surf.triangulate_hole()

    surf2 = Surface()
    surf2.make_cylinder(0., 0., -3., 0., 0., 3., 1.0, 21)
    maker = Domain([surf, surf2])

    maker.add_sharp_border_edges(surf)
    maker.set_borders()
    maker.create_mesh(9.0)
    maker.save("pvs_artery.mesh")
