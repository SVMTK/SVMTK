
import SVMTK as svm




if __name__ == "__main__":
   surf = svm.Surface() 

   surf.make_sphere(0.0,0.0,0.0,6.0,1.0)            


   slice_ = surf.slice(0,0,1.0,0.)                   
   slice_.simplify(0.9)   
 

   slice_.create_mesh(32)     
   slice_.save("test.vtu")
