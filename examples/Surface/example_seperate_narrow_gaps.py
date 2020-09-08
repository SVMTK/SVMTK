
import SVMTK as svm




if __name__ == "__main__":
    

   

   surf = svm.Surface("../Data/lh-pial.stl") 

   while(surf.separate_narrow_gaps()):
         surf.collapse_edges(0.3)
         #surf.smooth_taubin(1) 
   surf.save("seperate_close_junctures.off")
