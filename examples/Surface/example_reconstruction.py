""" Reconstruct a surface """

import SVMTK as  svm
from pathlib import Path

if __name__ == "__main__":
    print("Start ",__file__) 
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)

    surf = svm.Surface() 

    ### Reconstruction Parameters ###  
    angular_bound=20
    radius_bound=0.01
    distance_bound=0.01

    surf.reconstruct("../Data/kitten.xyz",angular_bound,radius_bound,distance_bound)


    surf.save(str(outdir/"reconstruct.stl"))
    print("Finish ",__file__)        
