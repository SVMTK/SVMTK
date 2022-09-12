"""  Separate overlapping surfaces """

import SVMTK as svm
from pathlib import Path
# TODO : improve 

if __name__ == "__main__":
   print("Start  ",__file__)
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)      
   sphere1 = svm.Surface() 
   sphere2 = svm.Surface() 
   
   # Create two spheres 
   # First sphere center coordinates, radius and edge length 
   x, y, z, r, edge_length = -1, 0, 0, 3.0, 0.2
   sphere1.make_sphere(x, y, z, r, edge_length) 
   
   # Second sphere center coordinates, radius and edge length  
   x, y, z, r, edge_length = 1, 0, 0, 3.0, 0.2
   sphere2.make_sphere(x, y, z, r, edge_length) 
   
   ### Separate overlapping surfaces parameters ### 
   # Vertex movement factor  
   edge_movment = -0.25     
   # Laplacian smoothing factor
   smoothing = 1.00
   # Number of iterations 
   max_iter = 100
   svm.separate_overlapping_surfaces(sphere1,sphere2,edge_movment,smoothing,max_iter) 
  
   ### Separate close surfaces parameters ### 
   # Vertex movement factor  
   edge_movment = -0.25    
   # Laplacian smoothing factor
   smoothing = 0.25
   # Number of iterations 
   max_iter = 100
   svm.separate_close_surfaces(sphere1,sphere2,edge_movment,smoothing,max_iter)  
   
   sphere1.save(str(outdir/"sphere1_disjointed.stl"))
   sphere2.save(str(outdir/"sphere2_disjointed.stl"))
   print("Finish ",__file__)   
      

