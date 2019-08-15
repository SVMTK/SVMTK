"""Convert a surface to a triangular mesh."""

import SVMTK as svm
if __name__ == "__main__":
    
 
   surf1 = svm.Surface()

   surf1.make_sphere(0,0,0,4) 

   p = svm.Point_3(10.0,0.0,0.0)
   surf1.extension(p,1.0,2.0,True)
   surf1.save("extension.stl")

