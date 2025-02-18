"""  Surface edits """

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start  ",__file__)
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)      
   sphere1 = svm.Surface() 
   sphere2 = svm.Surface() 
   sphere5 = svm.Surface() 
   # Create two spheres 
   # First sphere center coordinates, radius and edge length 
   x, y, z, r, edge_length = -1, 0, 0, 3.0, 0.2
   sphere1.make_sphere(x, y, z, r, edge_length) 
   sphere1.save(str(outdir/"init_embed1.stl"))   
   # Second sphere center coordinates, radius and edge length  
   x, y, z, r, edge_length = 1, 0, 0, 3.0, 0.2
   sphere2.make_sphere(x, y, z, r, edge_length) 
   sphere2.save(str(outdir/"init_embed2.stl"))
   
   sphere3 = svm.Surface(sphere1)
   sphere4 = svm.Surface(sphere1)
  
   x, y, z, r, edge_length = 1, 0, 0, 3.01, 0.2
   sphere5.make_sphere(x, y, z, r, edge_length) 
   sphere5.save(str(outdir/"init_embed3.stl"))

   
   ###  Embed  ### 
   # Embed function parameters 
   adjustment = -0.5 
   smoothing  =  0.2
   max_iter   = 100
   sphere1.embed(sphere2,adjustment,smoothing, max_iter) 
   sphere1.save(str(outdir/"embed.stl"))


   ###  Expose  ### 
   # Expose function parameters 
   adjustment = -0.5 
   smoothing  = 0.4
   max_iter   = 1000
   sphere3.expose(sphere2,adjustment,smoothing, max_iter)   
   sphere3.save(str(outdir/"expose.stl"))  
   

   ###  Enclose ### 
   # Enclose function parameters 
   adjustment =  0.5 
   smoothing  = 0.2
   max_iter   = 100
   sphere4.enclose(sphere2,adjustment,smoothing, max_iter)     
   sphere4.save(str(outdir/"enclose.stl")) 
 
        
   ### Separate ###
   # separate function parameters 
   adjustment = 0.5
   smoothing  = 0.4   
   max_iter   = 100
   sphere5.separate(sphere2,adjustment,smoothing, max_iter)     
   sphere5.save(str(outdir/"separate.stl"))

   print("Finish ",__file__)   
      

