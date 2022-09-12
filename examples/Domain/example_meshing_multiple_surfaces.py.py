""" Create a mesh with multiple surfaces """ 

import SVMTK as svm
import numpy as np
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)       
   surf = svm.Surface();
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   sphere1 = svm.Surface()
   sphere2 = svm.Surface()
   sphere3 = svm.Surface()    
     
   # Center coordinates    
   x, y, z = 0,0,0
   # Sphere radiuses 
   r1, r2, r3 = 5,10,15 
   # Edge length 
   edge_length = 1.
   # Create spheres 
   sphere1.make_sphere(x,y,z,r1,edge_length)
   sphere2.make_sphere(x,y,z,r2,edge_length)
   sphere3.make_sphere(x,y,z,r3,edge_length)

   maker = svm.Domain([sphere3,sphere2,sphere1])
   maker.create_mesh(22.)
   maker.save(str(outdir/"meshing_multiple_surfaces.mesh"))
   print("Finish ",__file__)       

