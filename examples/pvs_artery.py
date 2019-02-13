import brainmesh as bm


def pvs(x, y, z):
    if z**2 < 9:
        return x**2 + 4.*y**2 + x**2*y**2 - 4
    else:
        return z**2 - 9


if __name__ == "__main__":
    surf = bm.BrainSurface()
    surf.implicit_surface(pvs, 6.0, 30, 0.1, 0.1)

    surf2 = bm.BrainSurface()
    surf2.make_cylinder(0., 0., -3., 0., 0., 3., 1.0, 21)

    maker = bm.BrainMesh([surf, surf2])
    maker.add_sharp_border_edges(surf2)
    maker.set_borders()
    maker.default_creating_mesh()

    maker.save_mesh("pvs_artery.mesh")
