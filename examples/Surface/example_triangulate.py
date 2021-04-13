"""Convert a surface to a triangular mesh."""

import SVMTK as svm
from pathlib import Path
if __name__ == "__main__":
    print("Start ",__file__) 
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)

    p = svm.Surface("../Data/P.off")
    p.triangulate_faces()
    p.save(str(outdir/"triangularP.off"))
    print("Finish ",__file__)    
