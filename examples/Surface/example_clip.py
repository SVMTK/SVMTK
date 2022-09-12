"""Clips a triangular mesh."""

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
    print("Start ",__file__) 
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)

 
    surface = svm.Surface("../Data/lh-pial.stl")
    # The surface must be without self-intersections.
    surface.repair_self_intersections()
    # The surface should also not have holes.
    surface.fill_holes()  
    # Copy processed surface.
    surface2 = svm.Surface( surface )
     
        
    plane = svm.Plane_3(0,1,0,20) 
    preserve_manifold = True 
    surface.clip(plane, True)
    
    surface.save(str(outdir/"clipPlane.stl"))
    
    clipper = svm.Surface()
    clipper.make_cube(0,0,0,-80, -80 ,40, 2.0)
    
    # Parameters 
    invert = False 
    preserve_manifold = True 
    surface2.clip(clipper, invert, preserve_manifold)
    surface2.save(str(outdir/"clipSurf.stl"))  
    
      
    print("Finish ",__file__)    
