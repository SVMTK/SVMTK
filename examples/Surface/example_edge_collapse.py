""" Collapse long edges """

import SVMTK as svm
from pathlib import Path

if __name__ == "__main__":
   print("Start  ",__file__)
   outdir = Path("results")
   outdir.mkdir(exist_ok=True)
   
   shark = svm.Surface("../Data/mech-holes-shark.off")
   print("Number of edges before: {}".format( shark.num_edges()) )
   # Target edge length 
   target_edge_length = 0.2 
   shark.collapse_edges(target_edge_length)
   print("Number of edges after: {}".format(shark.num_edges()))
   shark.save(str(outdir/"collapsed_shark.stl"))
   
   # Create Cube 
   cube = svm.Surface()
   edge_length=0.2
   x0,y0,z0 = 0,0,0
   x1,y1,z1 = 4,4,4
   cube.make_cube(x0,y0,z0,x1,y1,z1,edge_length)
   print("Number of edges before: {}".format( cube.num_edges()) )   
   
   shark.collapse_edges()
   print("Number of edges before: {}".format( cube.num_edges()) )   
   surf.save(str(outdir/"collapsed_cube.stl"))

   print("Finish ",__file__)
