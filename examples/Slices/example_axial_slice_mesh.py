import SVMTK as svm
import numpy as np


if __name__ == "__main__":



   surf1 = svm.Surface()
   surf2 = svm.Surface()
    
   surf1.make_cube(0,0,0,5,5,5,0.5) 
   surf2.make_sphere(2.5,2.5,2.5,2.0,0.1)

   smap = svm.SubdomainMap()
   
   smap.add("10",2)
   smap.add("11",3)
   domain =svm.Domain([surf1,surf2],smap)
   domain.create_mesh(32)
   domain.remove_subdomain(3)
   

   surfaces = domain.get_boundaries()
   
   z = np.linspace(0.2,4.8,8) 




   for i in z:
     slice_ = svm.Slice(0,0,1,-i)
     slice_.slice_surfaces(surfaces)
     slice_.create_mesh(32.0)
     slice_.add_surface_domains([surf1,surf2],smap)
     slice_.remove_subdomain(3)
     slice_.save("axial-slice-{}.mesh".format(str(i) ) )


