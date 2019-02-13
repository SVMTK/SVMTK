import brainmesh as bm


if __name__ == "__main__":
    
 
   surf = bm.BrainSurface();

   surf.make_cube(-2.,-2.,-2.,2.,2.,2.)
  
   maker = bm.BrainMesh(surf)

   maker.add_sharp_border_edges(surf)

   line1 = [ bm.Point_3(0,0,-1.0),bm.Point_3(0,0,0.0), bm.Point_3(0,1.0,1.0)] 
   line2 = [ bm.Point_3(0,0,0.0), bm.Point_3(0,-1.0,1.0) ]
 
   maker.add_feature( line1)
   maker.add_feature( line2)

   maker.set_borders()
   maker.set_features()

   maker.create_mesh(24.)
   maker.excude(100, 0)
   maker.save_mesh("bifurcation_in_cube.mesh")
