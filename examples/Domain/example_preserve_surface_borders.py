import SVMTK as svm


if __name__ == "__main__":
   surf1 = svm.Surface()
   
   surf1.make_cube(0,0,0,20,20,20,1.5) 

   maker = svm.Domain(surf1)

   # Angle between edges that define sharp,
   angle = 70
   
   # Default angle is 60, at 90 some edges is fail 
   # to be detected  
   maker.add_sharp_border_edges(surf1, angle )
   maker.create_mesh(24) # need high resolution
   
   maker.save("box_edges.mesh")

