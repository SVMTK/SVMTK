
import SVMTK as svm


def plane(x,y,z) : 
    return z-3

if __name__ == "__main__":
    

   surf = svm.Surface();

   surf.implicit_surface(plane, 6.0,30,50,0.1)
   pial  =   svm.Surface("lh-pial.stl");

   pial.intersection(surf) 

   pial.save("test.stl")
   surf.save("plane.stl")


   
 


