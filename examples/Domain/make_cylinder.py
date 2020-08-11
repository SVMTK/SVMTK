
import SVMTK as svm


def pvs(x,y,z) : 
    return x**2 +4.*y**2 + x**2*y**2 -4.

if __name__ == "__main__":
    

   surf = svm.Surface();
   surf2 = svm.Surface();
 

   surf.make_cylinder(0.,0.,-3.,0.,0.,3.,0.8,60)
   surf2.make_cylinder(0.,0.,-3.,0.,0.,3.,0.7,60)

   surf2.intersection(surf)
   surf.difference(surf2) 
   surf.save("cylinder1.stl")
   surf2.save("cylinder2.stl")
   maker = svm.Domain(surf)
   
   maker.add_surface_points(surf)
   maker.add_surface_points(surf2)
   maker.refine_mesh(42)


   maker.save("cylinder.mesh")
