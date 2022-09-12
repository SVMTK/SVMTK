""" Creates a cone,cylinder,cube and sphere surface mesh """

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)    
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   surf = svm.Surface()
   
   # Edge length of surface mesh
   edge_length = 0.2

   ### Creating cube mesh ### 
   
   # Cube corner coordinates
   x0,y0,z0 = 0,0,0
   x1,y1,z1 = 4,4,4
   surf.make_cube(x0,y0,z0,x1,y1,z1,edge_length)
   surf.save(str(outdir/"cube.stl"))
   surf.clear()
   
   ### Creating sphere mesh ###
     
   # Sphere center coordinates
   x0,y0,z0 = 3,0,-2
   
   # Sphere radius 
   r = 2.5
   surf.make_sphere(x0,y0,z0,r,edge_length)
   surf.save(str(outdir/"sphere.stl"))
   surf.clear()   
 
   ### Creating cylinder mesh ###
   
   # Lower and upper cylinder center
   x0,y0,z0 = 0, 0, -3.0
   x1,y1,z1 = 0, 0, 3.0  
  
   # Cylinder radius 
   r = 1.0
   surf.make_cylinder(x0,y0,z0,x1,y1,z1,r,edge_length)
   surf.save(str(outdir/"cylinder.stl"))
   surf.clear()  
      
   ### Creating cone mesh ###
   
   # Lower and upper cone center  
   x0,y0,z0 = 0, 0, -3.0
   x1,y1,z1 = 0, 0, 3.0
   
   # Lower and upper radius
   r0,r1 = 3.0,0.2
   surf.make_cone(x0,y0,z0,x1,y1,z1,r0,r1,edge_length)
   surf.save(str(outdir/"cone.stl"))
   surf.clear()

   print("Finish ",__file__)
