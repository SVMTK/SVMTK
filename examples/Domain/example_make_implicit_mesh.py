""" Create a 3D mesh given an implicit function"""

import SVMTK as svm
from pathlib import Path

def chair_function ( x, y, z):
  x2=x*x 
  y2=y*y
  z2=z*z
  x4=x2*x2
  y4=y2*y2
  z4=z2*z2;
  return x4-1.2*x2*y2+3.6*x2*z2-7.50*x2+y4+3.6*y2*z2-7.50*y2+.2*z4-7.50*z2+64.0625-16.0*z*y2+16.0*x2*z

if __name__ == "__main__":
   print("Start ",__file__)   
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   surf = svm.Surface()
   # Implicit surface function parameters 
   bounding_sphere_radius = 6.0
   angular_bound = 30 
   radius_bound = 0.1
   distance_bound =0.1
   # Creates surface mesh given implicit function.
   surf.implicit_surface(chair_function, bounding_sphere_radius ,angular_bound ,radius_bound ,distance_bound )
   maker = svm.Domain(surf)
   maker.create_mesh(20)
   maker.save(str(outdir/"chair.mesh"))
   print("Start ",__file__)      
