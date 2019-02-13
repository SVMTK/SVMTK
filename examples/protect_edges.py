import brainmesh as bm


if __name__ == "__main__":
   surf = bm.BrainSurface();
   surf.make_cylinder(0., 0., 0., 2., 4., 3., 1.0, 21)

   #surf.split_edges(0.3)
   maker = bm.BrainMesh(surf)
   maker.add_sharp_border_edges(surf)
   maker.set_borders()
   maker.default_creating_mesh()
   #maker.refine_mesh()
   maker.save_mesh("cylinder.mesh")
