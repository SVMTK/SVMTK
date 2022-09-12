"""Compute boolean operations between blobby and eight."""

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
    print("Start  ",__file__)
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)

    surf1 = svm.Surface("../Data/blobby.off")
    surf2 = svm.Surface("../Data/eight.off")

    bop = svm.Surface(surf1)
    bop.difference(surf2)
    bop.save(str(outdir/"blobby_difference_eight.stl"))

    bop = svm.Surface(surf1)
    bop.union(surf2)
    bop.save(str(outdir/"blobby_union_eight.stl"))

    bop = svm.Surface(surf1)
    bop.intersection(surf2)
    bop.save(str(outdir/"blobby_intersection_eight.stl"))
    
    print("Finish ",__file__)
