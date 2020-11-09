
import SVMTK as svm


if __name__ == "__main__":
    

   surf = svm.Surface();
   surf2 = svm.Surface();
 

   surf.make_cylinder(0.,0.,-3.,0.,0.,3.,4,180)
   surf.save("cylinder.stl")
   maker = svm.Domain(surf)

   maker.add_sharp_border_edges(surf,60)

   maker.create_mesh(28)
   maker.exude(10,0) 

   maker.save("cylinder.mesh")
