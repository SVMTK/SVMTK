import brainmesh as bm


if __name__ == "__main__":
    sphere = bm.BrainSurface()
    sphere.make_cylinder(0., 0., 0., 2., 4., 3., 1.0, 21)

    sphere.split_edges(0.3)
    maker = bm.BrainMesh(sphere)
    maker.add_sharp_border_edges(sphere)
    maker.set_borders()
    maker.default_creating_mesh()
    # maker.refine_mesh()
    maker.save_mesh("results/cylinder.mesh")
