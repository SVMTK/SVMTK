""" Repairs self intersections """

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
    print("Start  ",__file__)   
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)
    pial = svm.Surface("../Data/lh-pial.stl")
    print("Number of self intersections", pial.num_self_intersections())

    ### Repair self-intersection parameters ### 
    
    # 
    volume_threshold = 0.01 
    
    #
    cap_threshold = 170
    
    #
    needle_threshold = 3.0
    
    #
    collapse_threshold = 0.2
    pial.repair_self_intersections(volume_threshold, cap_threshold, needle_threshold, collapse_threshold)
    print("Number of self intersections", pial.num_self_intersections() )

    pial.save(str(outdir/"repaired_pial.stl"))
    print("Finish ",__file__)
