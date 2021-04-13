"""Remesh a surface."""

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
    print("Start ",__file__)   
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)

    pig = svm.Surface("../Data/pig.off")
    target_edge_length = 0.04
    nb_iter = 3
    protect_border = True

    pig.isotropic_remeshing( target_edge_length,
                             nb_iter,
                             protect_border)

    pig.save(str(outdir/"remeshed_pig.off"))
    print("Finish ",__file__)
