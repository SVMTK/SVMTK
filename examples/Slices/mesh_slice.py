
import SVMTK as svm




if __name__ == "__main__":
    

   

   surf = svm.Surface() 

   surf.make_sphere(0.0,0.0,0.0,6.0)                 # Test surface 

   slice_ = surf.slice(0,0,1.0,0.)                   # Given that ax+by+cz+d then the inputs are a,b,c,d , returns a slice object
                                                     # The slice object contains X polylines which indicates the cutting of the surface.

   slice_.simplify(0.9)                              # simplifies all polylines, input is the stop criteria that (current count / initial count ) <= input
                                                     # i.e 0.5 removes half of polylines


   slice_.create_mesh(32)     
   slice_.save("test.off")
