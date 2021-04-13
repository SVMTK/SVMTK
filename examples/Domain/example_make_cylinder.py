
import SVMTK as svm


if __name__ == "__main__":
    

   surf = svm.Surface();
 

   surf.make_cylinder(0.,0.,-3.,0.,0.,3.,4,1.2)

   maker = svm.Domain(surf)

   maker.add_sharp_border_edges(surf,60)

   maker.create_mesh(28)
   maker.exude(10,0) 
   maker.boundary_segmentations()
   maker.save("cylinder.mesh")
