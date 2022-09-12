"""  """  

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)   
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)   
   surf = svm.Surface("../Data/lh-pial.stl") 

   # Create a cut plane 
   plane = svm.Plane_3(0,0,1,-5.41)
   # Cut surface 
   slce = surf.get_slice(plane)
   # Creating a mesh
   slce.create_mesh(32.0)
   # Saves mesh with helping faces
   slce.save("compact.stl")
   # Remove any helping faces
   slce.add_surface_domains([surf])
   slce.save("keep_all_components.vtu")    
   # Deletes smaller clusters of faces 
   slce.keep_largest_connected_component()
   slce.save("keep_largest_component.vtu")    

   print("Finish ",__file__)   
