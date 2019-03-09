
import SVMTK as svm




if __name__ == "__main__":
    

   

   surf = svm.Surface("lh-pial.stl") 

   surf.seperate_close_junctures(0.1)
   surf.save("seperate_close_junctures.off")
