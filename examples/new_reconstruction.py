import brainmesh as bm

eight = bm.BrainSurface("data/pig.off")
print(eight.self_intersections())
eight.reconstruct()
print(eight.self_intersections())

eight.save("foo.off")
