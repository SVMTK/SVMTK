
import SVMTK as svm




if __name__ == "__main__":
    

   surf = svm.Surface("lh-pial.stl");

   surf.isotropic_remeshing(1.0,2,False)
   surf.mesh_slice(0,0,1,0)

