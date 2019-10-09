import SVMTK as svm



if __name__ == "__main__":
    


   #surf.implicit_surface(chair_function, 6.0,30,0.1,0.1) #better
   surf = svm.Surface("../Data/lh-pial.stl");
       #surf.make_sphere(0.,0.,0.,1.0)
   #surf.isotropic_remeshing(1.0,2,False)

   for i in range(6):

       domain = svm.Domain(surf)
       domain.create_mesh(2**i)
       print(domain.number_of_cells())
  
