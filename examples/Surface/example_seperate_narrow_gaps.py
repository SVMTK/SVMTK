
import SVMTK as svm




if __name__ == "__main__":
   print("Start ",__file__)    

   surf = svm.Surface("../Data/lh-pial.stl") 
   print(surf.separate_narrow_gaps(-0.6,0.4,50)) 
   surf.save("seperate_close_junctures.off")
   print("Finish ",__file__)
