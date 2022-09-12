""" Create a four cube meshes  """

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)   
   surf = svm.Surface()
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   
   # Creating bounding cube 
   cube1 = svm.Surface()
   # Cube corner coordinates and edge length.
   x0,y0,z0,x1,y1,z1, edge_length = -4, -4, -4, 4, 4, 4, 0.5
   # Create a cube surface in 3D.
   cube1.make_cube(x0,y0,z0,x1,y1,z1,edge_length)
   cube1.save(str(outdir/"cube1.stl"))

   cube2 = svm.Surface()
   # Cube corner coordinates and edge length.
   x0,y0,z0,x1,y1,z1, edge_length = -2, -2, -2, 0, 2, 2, 0.5
   # Create a cube surface in 3D.
   cube2.make_cube(x0,y0,z0,x1,y1,z1,edge_length)
   cube2.save(str(outdir/"cube2.stl"))
   
   cube3 = svm.Surface()
   # Cube corner coordinates and edge length.
   x0,y0,z0,x1,y1,z1, edge_length = 0, -2, -2, 2, 2, 2, 0.5
   # Create a cube surface in 3D.
   cube3.make_cube(x0,y0,z0,x1,y1,z1,edge_length)      
   cube3.save(str(outdir/"cube3.stl"))   
   
   cube4 = svm.Surface()
   # Cube corner coordinates and edge length.
   x0,y0,z0,x1,y1,z1, edge_length = -1, -1, -1, 1, 1, 1, 0.5
   # Create a cube surface in 3D.
   cube4.make_cube(x0,y0,z0,x1,y1,z1,edge_length)         
   cube4.save(str(outdir/"cube4.stl"))   
   
   # Creating a SubdomainMap with 4 surface 
   smap = svm.SubdomainMap(4) 
   # Adding subdomains 
   smap.add("1000",1)
   smap.add("1010",2)
   smap.add("1100",3)   
   smap.add("*1",4)
   print(smap) 
 

   maker = svm.Domain([cube1,cube2,cube3,cube4],smap)
   # Angle between edges that define sharp,
   angle = 70
   # Default angle is 60, at 90 edges may fail 
   # to be detected  
   maker.add_sharp_border_edges(cube1, angle )
   maker.create_mesh(24) 
   maker.save(str(outdir/"four_cube.mesh"))   
   
   
