from svmtk import Surface

eight = Surface("data/pig.off")
print(eight.num_self_intersections())
eight.reconstruct()
print(eight.num_self_intersections())

eight.save("results/reconstructed_pig.off")
