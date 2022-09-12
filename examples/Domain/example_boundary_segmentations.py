""" Creates a cylinder mesh with boundary segmentations """
# TODO: complete 
import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start ",__file__) 
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   surf = svm.Surface()
   # Lower and upper cylinder center coordinates.
   x0,y0,z0,x1,y1,z1 = 0,0,-3,0,0,3
   # Cylinder radius and edge length. 
   radius, edge_length= 4, 1.2
   surf.make_cylinder(x0,y0,z0,x1,y1,z1,radius,edge_length)
   maker = svm.Domain(surf)
   maker.add_sharp_border_edges(surf,60)
   maker.create_mesh(28)
   maker.exude(10,0) 
   # The threshold angle to separate segmentation
   threshold_angle = 75 # default 85
   maker.boundary_segmentations(threshold_angle )
   maker.save(str(outdir/"cylinder.mesh")) 
   print("Finish ",__file__)    
