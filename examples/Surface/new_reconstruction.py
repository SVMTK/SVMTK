import SVMTK as svm

eight = svm.Surface("../Data/pig.off")
print(eight.self_intersections())
new_eight = svm.reconstruct_surface(eight)
print(new_eight.self_intersections())

new_eight.save("foo.off")
