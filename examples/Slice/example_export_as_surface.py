""" Export slice meshes as surface in 3D  """  

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)   
   outdir = Path("results")
   outdir.mkdir(exist_ok=True) 

   # Load Surface
   surf = svm.Surface("../Data/lh-pial.stl")
   
   ### Create and export slice mesh 1 ###
   
   # Plane equation parameters
   a,b,c,d = 0,0,1,-6.5
   # Construct a slice object with cut plane 
   slice_ = svm.Slice(a,b,c,d)
   # Cut multiple surfaces
   slice_.slice_surfaces([surf])
   slice_.create_mesh(32.0)
   slice_.add_surface_domains([surf])
   # Export slice as surface 
   surf2 = slice_.export_as_surface()
   surf2.save("slice_1.stl")

   ### Create and export slice mesh 2 ###

   # Plane equation parameters
   a,b,c,d = 0,0,1,-30.5
   # Construct a slice object with cut plane 
   slice_ = svm.Slice(a,b,c,d)
   # Cut multiple surfaces
   slice_.slice_surfaces([surf])
   slice_.create_mesh(32.0)
   slice_.add_surface_domains([surf])
   # Export slice as surface 
   surf2 = slice_.export_as_surface()
   surf2.save("slice_2.stl")
   
   print("Finish ",__file__)   
