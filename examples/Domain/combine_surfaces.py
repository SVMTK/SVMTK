
import SVMTK as svm


def pvs(x,y,z) : 
    return x**2 +4.*y**2 + x**2*y**2 -4.

if __name__ == "__main__":
    

   surf = svm.Surface();

   surf.implicit_surface(pvs, 6.0,30,0.1,0.1)
   
   surf.clip(0,0,-1,-3, False)
   surf.clip(0,0,1,-3, False)
   surf.triangulate_hole()

   surf2 = svm.Surface();

   surf2.make_cylinder(0.,0.,-3.,0.,0.,3.,0.9,60)


   maker = svm.Domain([surf,surf2])

   maker.add_sharp_border_edges(surf2,90)
   maker.add_sharp_border_edges(surf,90)
   maker.set_borders()

   maker.create_mesh(126)


   maker.save("pvs_artery.mesh")


   
 



