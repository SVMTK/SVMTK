""" Slice a surface in different planes  """  

import SVMTK as svm
from pathlib import Path
import numpy as np

if __name__ == "__main__":
   print("Start ",__file__)   
   cube  = svm.Surface()
   sphere= svm.Surface()
   outdir = Path("results")
   outdir.mkdir(exist_ok=True) 
      
   surf = svm.Surface("../Data/boundary.stl")
    
   # Get the span of the surface in z-direction.
   zmin, zmax = surf.span(2) 
   z = np.linspace(zmin,zmax,8) 

   for i in z[2:-2]:
     slce = boundary.slice(0,0,1,-i)
     slce.create_mesh(32.0)
     slce.save(f"slice-{i}.vtu")    

   print("Finish ",__file__)   
