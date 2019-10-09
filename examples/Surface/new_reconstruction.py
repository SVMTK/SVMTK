import SVMTK as svm

eight = svm.Surface("../Data/pig.off")
print(eight.num_self_intersections())
eight.reconstruct(10.,8.0,0.25,0.02,5.0)
print(eight.num_self_intersections())
eight.fill_holes()
print(eight.num_self_intersections())
eight.save("foo.off")
