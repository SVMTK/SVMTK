""" Combine slices from different surfaces  """  

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)   
   outdir = Path("results")
   outdir.mkdir(exist_ok=True) 

   cube  = svm.Surface()
   sphere= svm.Surface()

   # Cube corner coordinates 
   x0,y0,z0,x1,y1,z1 =0,0,0,5,5,5
   # Sphere center and radius  
   cx,cy,cz,r = 2.5, 2.5, 2.5,3.0
   # Edge length of the surfacemesh
   edge_length = 0.2
   
   # Create cube surface.    
   cube.make_cube(x0,y0,z0,x1,y1,z1,edge_length)
   # Create sphere surface. 
   sphere.make_sphere(cx,cy,cz,r,edge_length)

   # Point on plane and plane vector  
   point, vector= svm.Point_3(0,0,1.6),svm.Vector_3(0,0,1)
   plane = svm.Plane_3(point,vector)
   
   # Construct a slice object with cut plane 
   slice1 = cube.get_slice(plane)
   slice2 = sphere.get_slice(plane)
   slice1.add_constraints(slice2)

   slice1.create_mesh(8)
   slice1.save("combined_slices.vtu")

   print("Finish ",__file__)   
