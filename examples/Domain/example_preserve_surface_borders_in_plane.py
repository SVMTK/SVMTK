""" Create a cube mesh with perserved sharp edges in a plane """
# TODO USE POLYLINES 
import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__)   
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   surf = svm.Surface()
   # Center coordinates 
   x,y,z = 0,0,0 
   # Sphere radius 
   r = 5 
   #Edge length 
   edge_length = 0.2 
   surf.make_sphere(x,y,z,r,edge_length) 
   # Define a plane with a point on the plane and the normal vector.
   plane = svm.Plane_3(svm.Point_3(0,0,2), svm.Vector_3(0,0,1)) 
   surf.clip(plane)
   #surf.repair_self_intersections()
   surf.collapse_edges()
   surf.isotropic_remeshing(0.15,20,True)
   surf.save(str(outdir/"preserved_edge_in_plane.stl"))
   # Clips the surface with a circle in the plane with a given radius 
   radius = 8
   #surf.clip(svm.Point_3(0,0,2), svm.Vector_3(0,0,1),radius)

   maker = svm.Domain(surf)
   # Default angle is 60, at 90 edges may fail 
   # to be detected  
   angle = 40
   maker.add_sharp_border_edges(surf, plane, angle)
   maker.create_mesh()  
   maker.save(str(outdir/"preserved_edge_in_plane.mesh"))
   print("Finish ",__file__)   
