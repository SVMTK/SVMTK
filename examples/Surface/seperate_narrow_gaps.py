
import SVMTK as svm




if __name__ == "__main__":
    

   

   surf = svm.Surface("../Data/lh-pial.stl") 

   surf.separate_narrow_gaps()
   surf.save("seperate_close_junctures.off")
