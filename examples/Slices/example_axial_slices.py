import SVMTK as svm
import numpy as np


if __name__ == "__main__":
    


   boundary = svm.Surface("../Data/boundary.stl")
    
   zmin, zmax = boundary.span(2) 
   z = np.linspace(zmin,zmax,8) 

   for i in z[2:-2]:
     rhslice = boundary.slice(0,0,1,-i)
     rhslice.create_mesh(32.0)
     rhslice.save("axial-slice-{}.vtu".format(str(i) ) )


