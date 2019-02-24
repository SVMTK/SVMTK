
import SVMTK as svm


def pvs(x,y,z) : 


       return x**2 +4.*y**2 + x**2*y**2 -4.

if __name__ == "__main__":
    

   surf = svm.Surface();

   surf.implicit_surface(pvs, 6.0,30,0.1,0.1)
   
   surf.fill_holes()

   surf.save("fill_holes.off")
 



