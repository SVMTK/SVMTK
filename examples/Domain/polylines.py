"""Compute number of self-intersections."""

import SVMTK as  svm
from numpy import random as rand 
from pathlib import Path

outdir = Path("../results")
outdir.mkdir(exist_ok=True)

sphere = svm.Surface("../Data/sphere_r6.off") 

mesh = svm.Domain(sphere)
print(rand.random())
dx = 0.75
vector1 = [ svm.Point_3(rand.random(),rand.random(),i*dx-2.5) for i in range(6) ] 


mesh.add_feature(vector1)
mesh.set_features()
mesh.create_mesh(1)
mesh.save("polylinewithsphere.mesh")

