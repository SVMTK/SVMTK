""" Create a 3D mesh with conforming bifurcation """

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)   
   surf = svm.Surface()
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   
   # Cube corner coordinates and edge length.
   x0,y0,z0,x1,y1,z1, edge_length = -2, -2, -2, 2, 2, 2, 0.5
   # Create a cube surface in 3D.
   surf.make_cube(x0,y0,z0,x1,y1,z1,edge_length)
   
   maker = svm.Domain(surf)  
   
   # Edges that Angles to detect 
   angle_detect = 85
   maker.add_sharp_border_edges(surf,85)

   # Create 3 lines with the bifurcation point. 
   line0 = [svm.Point_3(0,0,0), svm.Point_3(0,1.0,1.0)] 
   line1 = [svm.Point_3(0,0,-1.0), svm.Point_3(0,0,0)] 
   line2 = [svm.Point_3(0,0,0), svm.Point_3(0,-1.0,1.0)]
   
   # Add features to the mesh construction.
   maker.add_feature(line0)
   maker.add_feature(line1)
   maker.add_feature(line2)

   # Create the mesh 
   maker.create_mesh(16.) 
   maker.exude(100, 0)
   print( maker.get_features( ) )
   maker.save(str(outdir/"bifurcation_in_cube.mesh"))
   print("Finish ",__file__)   
   
