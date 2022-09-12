""" Separate narrow gaps in the surface mesh """ 
# TODO: complete 
import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)    
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   surf = svm.Surface("../Data/lh-pial.stl")
   
   print( surf.num_self_intersections())
   
   #surf.repair_self_intersections()
   #surf.isotropic_remeshing(1.0,5,False)
   surf.save(str(outdir/"separate_narrow_gaps_init.stl"))
   ### Parameter inputs ###      
   # Multiplier of edge length that determines vertex displacment.
   adjustment = -1.0
   # Laplacian smoothing factor after each iteration
   smoothing = 0.0
   # Number of maximum iterations.
   max_iter = 100 
   
   surf.separate_narrow_gaps(adjustment,smoothing,max_iter) 
   print( surf.num_self_intersections() )
   surf.save(str(outdir/"separate_narrow_gaps.stl"))
   
   print("Finish ",__file__)
