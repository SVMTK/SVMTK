
import SVMTK as svm




if __name__ == "__main__":
    

   surf = svm.Surface()
   surf.make_sphere(0.,0.,0.,4)

   surf.mesh_slice(0,0,1,0,"mesh_slice.off")

