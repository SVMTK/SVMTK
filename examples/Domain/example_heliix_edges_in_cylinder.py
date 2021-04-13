import SVMTK as svm
import numpy as np

if __name__ == "__main__":
    


   lines =[]
    
   surf = svm.Surface();

   surf.make_cylinder(0.,0.,0.,0.,0.,10.,5.,1.0)
   R=3.0
   helix1=[]
   helix2=[]
   for i in np.linspace(0.5,9.5,40) : 
       helix1.append(svm.Point_3(   R*np.cos( np.pi*i*0.5), R*np.sin( np.pi*i*0.5) ,i) )
       helix2.append(svm.Point_3(  -R*np.cos( np.pi*i*0.5), -R*np.sin( np.pi*i*0.5) ,i)  ) 

        
   
   lines.append(helix1)
   lines.append(helix2)

   maker = svm.Domain(surf)


   for k in lines:
       maker.add_feature(k)
   
   maker.add_sharp_border_edges(surf,70)
   maker.create_mesh(22.)
   print(maker.get_curves())
   maker.boundary_segmentations()
   maker.save("Helixfeature.mesh")


