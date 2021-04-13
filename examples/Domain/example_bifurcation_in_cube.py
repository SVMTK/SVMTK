import SVMTK as svm


if __name__ == "__main__":
    
 
   surf = svm.Surface();

   surf.make_cube(-2.,-2.,-2.,2.,2.,2.,0.5)
   surf.save("cube_.off")
   maker = svm.Domain(surf)
  
   maker.add_sharp_border_edges(surf,85)

   line0 = [svm.Point_3(0,0,0.0), svm.Point_3(0,1.0,1.0)] 
   line1 = [ svm.Point_3(0,0,-1.0),svm.Point_3(0,0,0.0)] 
   line2 = [ svm.Point_3(0,0,0.0), svm.Point_3(0,-1.0,1.0) ]
 
   maker.add_feature( line0)
   maker.add_feature( line1)
   maker.add_feature( line2)

   maker.create_mesh(24.)
   maker.exude(100, 0)
   maker.save("bifurcation_in_cube.mesh", True)


   
