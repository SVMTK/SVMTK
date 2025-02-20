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

   smap = svm.SubdomainMap()
   # Adding subdomains 
   smap.add("100",1)
   smap.add("110",2)   
   smap.add("111",3)
   

   maker = svm.Domain([sphere3,sphere2,sphere1],smap)
   maker.create_mesh(32)
   
   bo1 = maker.get_interface( (1,0)) 
   bo2 = maker.get_boundary(1)
   bo1.save("bo1.stl")
   bo2.save("bo2.stl")
   print(maker.number_of_facets()  ) 
   import time
   start_time = time.time()
   maker.set_collision_spheres(1,0)
   print("--- %s seconds ---" % (time.time() - start_time))

   maker.remove_subdomain(1)
   maker.write_facet_data("facetdata_2.dat")
   print(maker.number_of_facets() )
   
   """   
   We want the distance from boundary to another surface.
   
   TODO : 
        
   Abstract class -> call 
   
   
   
   Solution 
   vertex attribute = double "depth" value  
   
   Problems : 
   
   
   closest surface in normal direction : 
   
   pial 1 
   pial 2 
   bounding 
   
   
   Everthing 
   
   
   
   
   
   """
   
   
   
   
   
   maker.save("reduced.mesh")
   print("Finish ",__file__)       

