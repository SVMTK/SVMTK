"""Smooth a surface."""

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
    print("Start ",__file__) 
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)

    pig1 = svm.Surface("../Data/pig.off")
    pig1.fill_holes()
    
    # Creating a deep copy
    pig2 = svm.Surface(pig1)
    pig3 = svm.Surface(pig1)
    
    # Laplacian smoothing factor and number of iterations : 
    smoothing_factor, num_iter = 0.8, 5
    pig1.smooth_laplacian(smoothing_factor, num_iter)

    # Taubin smoothing iterations
    num_iter = 5 
    pig2.smooth_taubin(num_iter)
    
    # Smooth shape speed (1e-6,1] and number of iterations
    speed, num_iter = 1.e-4,5
    pig3.smooth_shape(speed,num_iter)
    
    pig1.save(str(outdir/"Laplacian_smooth.stl"))
    pig2.save(str(outdir/"Taubin_smooth.stl")) 
    pig3.save(str(outdir/"Shape_smooth.stl")) 
    print("Finish ",__file__) 
