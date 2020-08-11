
import SVMTK as svm

import numpy as np


if __name__ == "__main__":
    

   #rhwhite = svm.Surface("rh-white.stl")



   #smap = svm.SubdomainMap()

   #smap.add("100",1)
   #smap.add("110",1)

   #domain = svm.Domain([rhpial,rhwhite,rhvent], smap) 


   #domain.create_mesh(16) 
   #boundary = domain.get_boundary()
   #boundary.save("bundary.stl")


   boundary = svm.Surface("../Data/boundary.stl")
    
   zmin, zmax = boundary.span(2) # should be in c++ code
   print(zmin,zmax)
   z = np.linspace(zmin,zmax,8) 

   for i in z[2:-2]:
     rhslice = boundary.slice(0,0,1,-i)
     rhslice.simplify(0.8)
     rhslice.create_mesh(32.0)
     rhslice.save("axial-slice-{}.vtu".format(str(i) ) )

       
   #surf2 = svm.Surface("test-output-slice.off")

   #surf2.isotropic_remeshing(1.8,2,False) 

   #surf2.save("afgf.off")
