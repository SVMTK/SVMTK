
import SVMTK as svm




if __name__ == "__main__":
   print("Start ",__file__)
   surf1 = svm.Surface() 
   surf2 = svm.Surface() 

   surf1.make_sphere(-1.0,0.0,0.0,3.0,0.2) 
   surf2.make_sphere(1.0,0.0,0.0,3.0,0.2) 
   surf1.save("surface_overlapp_before1.off")
   surf2.save("surface_overlapp_before2.off")
   

   print( svm.separate_overlapping_surfaces(surf1,surf2,-0.5,0.5,800) )

   print( svm.separate_close_surfaces(surf1,surf2,-0.5,0.8,100)  )
   print("Finish ",__file__)   
      

