"""Convert a surface to a triangular mesh."""

import SVMTK as svm
if __name__ == "__main__":
    print("Start ",__file__)   
 
   surf1 = svm.Surface()

   surf1.make_sphere(0,0,0,4,.4) 

   p = svm.Point_3(10.0,0.0,0.0)
   surf1.extension(p,0.3,8.0, 0.1,True)
   surf1.save("extension.stl")
   print("Finish ",__file__)

