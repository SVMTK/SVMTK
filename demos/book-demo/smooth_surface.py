#!/usr/bin/python

if __name__ =='__main__':
        import SVMTK as svm
        import argparse
        import os
        hemi , surface = os.path.basename(__file__).split(".")[0:2]
        parser = argparse.ArgumentParser(description='The smoothing of an input surface using Laplacian or volume preserving Laplacian'.format(hemi,surface) )


        parser.add_argument("--i"   ,type=str, help = "The path to the input surface.",required=True)
        parser.add_argument("--iter"  ,type=int , default= 1 ,help = "The number of smoothing iterations")
        parser.add_argument("--mdisp"  ,type=float , default= 1.0 ,help = "Scalar factor of the displacement vector of a point.")
        parser.add_argument("--volume" ,action='store_true',help = "Include to use volume preserving Laplacian smoothing, known as taubian smoothing")
        parser.add_argument("--o"   ,type=str, help = "The path to the output surface.")
        Z = parser.parse_args() 

        surf  = svm.Surface(Z.i)

        if not Z.volume : 
             surf.smooth_laplacian(Z.mdisp,Z.iter)
        else : 
             surf.smooth_taubin(Z.iter)


        surf.save(Z.o) # make wrtier





