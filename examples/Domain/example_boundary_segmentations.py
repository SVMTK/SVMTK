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
   x0,y0,z0,x1,y1,z1 = 0,0,-3,0,0,8
   # Cylinder radius and edge length. 
   radius, edge_length= 5, 0.707
   surf.make_cylinder(x0,y0,z0,x1,y1,z1,radius,edge_length)
   #surf.isotropic_remeshing(1.0,5,True)
   
   surf.save("tst.stl")

   maker = svm.Domain(surf)
   maker.add_sharp_border_edges(surf,60)
   
   es = .707
   fs = 0.707
   cs = .707
   fa = 30. 
   fd = 0.07
   r  = 4
       
   #maker.create_mesh(edge_size=es,cell_size=cs,facet_size=fs,facet_angle=fa,facet_distance=fd,cell_radius_edge_ratio=r) 
   
   maker.create_mesh(32.0)
   #maker.perturb()
   #maker.odt() 
   #maker.exude() 
   # The threshold angle to separate segmentation
   threshold_angle = 60 # default 85
   #print( maker.get_subdomains())
   print(maker.get_patches())   
   #maker.boundary_segmentations(1,threshold_angle)

   maker.save(str(outdir/"cylinder.mesh")) 
   print("Finish ",__file__)    
