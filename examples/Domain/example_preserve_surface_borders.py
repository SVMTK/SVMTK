""" Create a cube mesh with perserved sharp edges """

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)   
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   surf1 = svm.Surface()  
   # Creates a cube surfaces by specific two corners of the cube 
   # Corner coordinates 
   x0,y0,z0 = 0,0,0
   x1,y1,z1 = 20,20,20
   # Edge length
   edge_length = 1.5
   surf1.make_cube(x0,y0,z0,x1,y1,z1,edge_length) 
   maker = svm.Domain(surf1)
   # Angle between edges that define sharp,
   angle = 70
   # Default angle is 60, at 90 edges may fail 
   # to be detected  
   maker.add_sharp_border_edges(surf1, angle )
   # Need for relative high resolution
   maker.create_mesh(24) 
   maker.save(str(outdir/"box_edges.mesh"))
   print("Start ",__file__)   
