#!/usr/bin/python

if __name__ =='__main__':
        import SVMTK as svm
        import argparse
        import os

        parser = argparse.ArgumentParser(description='The smoothing of an input surface using Laplacian or volume preserving Laplacian')

        parser.add_argument("--i"   ,type=str, help = "The path to the input surface.",required=True)
        parser.add_argument("--iter"  ,type=int , default= 2 ,help = "The number of smoothing iterations. The default value is 2 .")
        parser.add_argument("--mdisp"  ,type=float , default= 1.0 ,help = "Scalar factor of the displacement vector of a point. The value should be <1.0, and the default value is 1.0")
        parser.add_argument("--volume" ,action='store_true',help = "Include to use volume preserving Laplacian smoothing, known as taubian smoothing")
        parser.add_argument("--o"   ,type=str, help = "The path to the output surface.")
        Z = parser.parse_args() 

        surf  = svm.Surface(Z.i)                       # Loads the input surface into the Surface Obejct

        if not Z.volume :                              
             surf.smooth_laplacian(Z.mdisp,Z.iter)     # Laplacian smoothing on the input surface 
        else : 
             surf.smooth_taubin(Z.iter)                # Taubin smoothing on the input surface 


        surf.save(Z.o)                                 # Saves the smooth surface.





