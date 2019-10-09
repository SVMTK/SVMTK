
import SVMTK as svm


def pvs(x,y,z) : 
    return x**2 +4.*y**2 + x**2*y**2 -4.

if __name__ == "__main__":
    

   surf = svm.Surface();

 

   surf.make_cylinder(0.,0.,-3.,0.,0.,3.,0.8,60)

   surf.save("cylinder.stl")
   maker = svm.Domain(surf)



   maker.add_corners(surf)

   maker.create_mesh(20)


   maker.save("cylinder.mesh")
