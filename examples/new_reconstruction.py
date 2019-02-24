import SVMTK as svm

eight = svm.Surface("data/pig.off")
print(eight.self_intersections())
new_eight = svm.reconstruct_surface(eight)
print(new_eight.self_intersections())

new_eight.save("foo.off")
