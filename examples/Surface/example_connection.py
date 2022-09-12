"""Convert a surface to a triangular mesh."""

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start  ",__file__)
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)     
   
   surf2 = svm.Surface()
   surf1 = svm.Surface()
      
   # Surface mesh edge length:   
   edge_length = 0.4    
   
   # First surface sphere with center coordinates and radius: 
   x0, y0, z0, r0 = 0, 0, 8, 4   
   surf2.make_sphere(x0,y0,z0,r0,edge_length) 
   
   # Second surface sphere with center coordinates and radius: 
   x1, y1, z1, r1 = 0, 0, -8, 4  
   surf1.make_sphere(x1,y1,z1,r1,edge_length) 

   # Creating a cylinder connection with radius between surfaces.
   cylinder_radius = 2.0
   bridge = surf1.connection(surf2, cylinder_radius, edge_length)
   
   # Combining the surfaces into a single surface mesh.
   bridge.union(surf1)
   bridge.union(surf2)

   bridge.save(str(outdir/"connected.stl"))
   print("Finish ",__file__)
