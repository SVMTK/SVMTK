""" Creates a cone,cylinder,cube and sphere surface mesh """

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)    
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   surf = svm.Surface()
   ### Defining an implicit surface ### 
   def torus(x,y,z):
       import numpy as np
       return (x**2 + y**2 +z**2 + 5**2 -  2**2)**2 - 4*5**2*(x**2+y**2)
 
   # Implicit surface fucntion parameters  
   bounding_sphere_radius = 6.0
   angular_bound = 30 
   radius_bound = 0.1
   distance_bound =0.1
   # Creates surface mesh given implicit function.
   surf.implicit_surface(torus,bounding_sphere_radius ,angular_bound ,radius_bound ,distance_bound ) 
   surf.save(str(outdir/"torus.stl"))
   print("Finish ",__file__)
