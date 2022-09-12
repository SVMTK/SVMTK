"""Convert a surface to a triangular mesh."""

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start  ",__file__)   
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)   
   sphere = svm.Surface()
   
   # Surface mesh edge length:
   edge_length = 0.2
   
   # Sphere center coordinates and radius:
   x,y,z,r =0,0,0,4
   sphere.make_sphere(x,y,z,r,edge_length) 
  
   # Point coordinates:
   px,py,pz = 10,0,0 
   p = svm.Point_3(px,py,pz)
  
   # Extension radius,length and edge length: 
   radius, length, ext_edge_length = 0.6,8.0,0.1
  
   # Enforce extension to be normal to the surface
   use_normal = True 
   extension = sphere.extension(p, radius, length, edge_length,use_normal)
   extension.save("temp.stl")
   extension.union(sphere)
   extension.save(str(outdir/"sphere_with_extension.stl"))
   print("Finish ",__file__)

