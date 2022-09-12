""" Create 3D mesh with two conforming edge helixes inside """ 

import SVMTK as svm
import numpy as np
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)       
   surf = svm.Surface();
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   
   cube = svm.Surface()
   sphere = svm.Surface()
   cylinder = svm.Surface()    
     
   # Edge length.
   edge_length = 1.     
     
   ### Create cube ###
   # Cube corner coordinates.
   x1, y1, z1, x2, y2, z2 = -1.5,-1.5,-1.5,1.5,1.5,1.5
   cube.make_cube(x1,y1,z1,x2,y2,z2,edge_length)

   ### Create sphere ###
   # Sphere center coordinates and radius. 
   x, y, z, r = 0.,0.,0.,5.
   sphere.make_sphere(x,y,z,r,edge_length)

   ### Create cylinder ###
   # Cylinder radius, top and bottom center coordinates. 
   x1, y1, z1, x2, y2, z2, r = -3,0,0,3,0,0,3
   cylinder.make_cylinder(x1,y1,z1,x2,y2,z2,r,edge_length)

   # Meshing without subdomainmap 

   # Create subdomainmap
   smap = svm.SubdomainMap()
   # Adding subdomains 
   smap.add("100",1)
   smap.add("110",2)
   smap.add("111",3)
   
   # List of surfaces, cylinder in cube, and cube in sphere. 
   surfaces = [sphere, cylinder, cube]
   # Meshing with subdomainmap.
   maker = svm.Domain(surfaces,smap)
   maker.add_sharp_border_edges(cylinder, 60 )
   maker.add_sharp_border_edges(cube, 60 )      
   maker.create_mesh(22.)

   maker.save(str(outdir/"mesh_without_interface_marking.mesh"))
   print( maker.get_patches() )
   # The interface tags will change how the mesh is written to file.
   smap.add_interface((1,0),8)
   smap.add_interface((2,3),16)
   smap.add_interface((1,2),32)    
   
   maker.save(str(outdir/"mesh_with_interface_marking.mesh"))
   print("Finish ",__file__)       

