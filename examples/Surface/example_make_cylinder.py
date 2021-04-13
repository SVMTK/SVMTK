
import SVMTK as svm


if __name__ == "__main__":
   print("Start ",__file__)    

   surf = svm.Surface();
   surf2 = svm.Surface();
 

   surf.make_cylinder(0.,0.,-3.,0.,0.,3.,0.2,0.1)
   surf.save("cylinder.stl")


   surf2.make_cube(-1,-1,-1,1,1,1,0.2)
   surf2.save("cube.stl")   
   surf2.clip(surf, invert=True)
   surf2.save("re.stl")
   print("Finish ",__file__)
