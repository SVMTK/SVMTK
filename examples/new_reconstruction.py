import brainmesh as bm

eight = bm.BrainSurface("data/pig.off")
print(eight.self_intersections())
new_eight = bm.reconstruct_surface(eight)
print(new_eight.self_intersections())

new_eight.save("foo.off")
