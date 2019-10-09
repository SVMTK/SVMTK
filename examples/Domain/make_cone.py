
import SVMTK as svm


if __name__ == "__main__":
    

   surf = svm.Surface();

   surf.make_cone(0.,0.,-3.,0.,0.,3.,0.9, 0.0, 60)

   surf.save("cone.stl")



   
 



