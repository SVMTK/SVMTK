""" Add subdomain tags to a 2D mesh """  

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)   
   outdir = Path("results")
   outdir.mkdir(exist_ok=True) 
   cube  = svm.Surface()
   sphere1 = svm.Surface()
   sphere2 = svm.Surface() 
   # Edge length of surface mesh
   edge_length = 0.2      
   # Cube corner coordinates 
   x0,y0,z0,x1,y1,z1 =0,0,0,5,5,5
   # Sphere center coordinates and inner radius and outer radius 
   cx,cy,cz,r1,r2 = 2.5, 2.5, 2.5, 2.5,3.54

   cube.make_cube(x0,y0,z0,x1,y1,z1,edge_length)
   sphere1.make_sphere(cx,cy,cz,r1,edge_length)
   sphere2.make_sphere(cx,cy,cz,r2,edge_length)
   
   smap = svm.SubdomainMap() 
   # Sphere2 enclose cube, and cube encloses sphere1
   surfaces = [sphere2,cube,sphere1]
   # Marks cells only inside the outer sphere as 1 
   smap.add("100",1)
   # Marks cells inside cube and outer sphere as 2 
   smap.add("110",2)
   # Marks cells inside all surfaces as 3 
   smap.add("111",3)
   # Plane equation parameters
   a,b,c,d = 0,0,1,-2.5
   # Construct a slice object with cut plane 
   slice_ = svm.Slice(a,b,c,d)
   # Cut multiple surfaces
   slice_.slice_surfaces(surfaces)
   slice_.create_mesh(32.0)
   slice_.add_surface_domains(surfaces,smap)
   slice_.save(str(outdir/"slice_with_3subdomains.mesh") )
   slice_.remove_subdomain(3)
   slice_.save(str(outdir/"slice_with_2subdomains.mesh") )
   
   print("Finish ",__file__)   
