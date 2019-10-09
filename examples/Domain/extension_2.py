"""Convert a surface to a triangular mesh."""

import SVMTK as svm
if __name__ == "__main__":
    
 
   surf1 = svm.Surface("../Data/ernie-ventricles.stl")
   surf1.smooth_laplacian(0.6,3) 
   surf1.isotropic_remeshing(0.3,2,False)
   surf1.smooth_taubin(4) 
   surf1.fill_holes()
   p = svm.Point_3(0.00,-37.33,-61.53)
   surf1.extension(p,1.0,10.0,True)
   
  
   surf1.save("extension_2.stl")
   

   domain =svm.Domain(surf1) 
   #domain.add_corners(surf1)
   domain.add_sharp_border_edges(surf1,80) 
   domain.set_borders()


   domain.create_mesh(168)

   domain.save("extension3.mesh")
