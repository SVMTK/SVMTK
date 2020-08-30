
import SVMTK as svm


def pvs(x,y,z) : 
    return x**2 +4.*y**2 + x**2*y**2 -4.

if __name__ == "__main__":
    

   surf = svm.Surface();
   surf2 = svm.Surface();
 

   surf.make_cylinder(0.,0.,-3.,0.,0.,3.,0.8,60)

   maker = svm.Domain(surf)

   maker.add_sharp_border_edges(surf,90)

   maker.create_mesh(42)


   maker.save("cylinder.mesh")
