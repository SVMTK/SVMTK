import unittest
import SVMTK
def ellipsoid_function( x, y, z):
  return x*x + 4.*y*y +4.*z*z-1.;

class Surface_Test(unittest.TestCase):



    def test_surface_io(self):
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,1.,1.,1.,1) 
        surface.save("cube.off")
        surface.save("cube.stl")
        surface= SVMTK.Surface("cube.off")
        self.assertTrue(surface.num_vertices()==8 and surface.num_faces()==12 and surface.num_edges()==18) 
        surface= SVMTK.Surface("cube.stl") 
        self.assertTrue(surface.num_vertices()==8 and surface.num_faces()==12 and surface.num_edges()==18)


    def test_surfaces_shapes(self):
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,1.,1.,1.,1) 
        self.assertTrue(surface.num_vertices()==8 and surface.num_faces()==12 and surface.num_edges()==18) 
        surface.clear() 
        self.assertEqual( surface.num_vertices(), 0) 
        surface.make_cone(0.,0.,0.,0,0.,2., 4.0,2.0,3)
        self.assertTrue(surface.num_vertices()==8 and surface.num_faces()==12 and surface.num_edges()==18) 
        surface.clear() 
        self.assertEqual( surface.num_vertices(), 0) 
        surface.make_cone(0.,0.,0.,0,0.,2., 1.0,0.0,3) 
        self.assertTrue(surface.num_vertices()==20 and surface.num_faces()==36 and surface.num_edges()==54) 
        surface.clear() 
        self.assertEqual( surface.num_vertices(), 0) 
        surface.make_cylinder(0.,0.,0.,1.,1.0,.1,2.0,4) 
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36) 
        surface.clear() 
        self.assertEqual( surface.num_vertices(), 0) 




    def test_surface_remeshing(self):
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,1.,1.,1.,1) 
        self.assertTrue(surface.num_vertices()==8 and surface.num_faces()==12 and surface.num_edges()==18) 
        surface.isotropic_remeshing(1.0,1,1)
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36)

    def boolean_operations(self):
        surface1 =SVMTK.Surface()  
        surface1.make_cube(0.,0.,0.,1.,1.,1.,1) 
        surface2 =SVMTK.Surface()  
        surface2.make_cube(0.5,0.5,0.5,1.5,1.5,1.5,1) 
        surface1.union(surface2) 
        self.assertTrue(surface.num_vertices()==20 and surface.num_faces()==36 and surface.num_edges()==54) 
        surface1.clear() 
        surface1.make_cube(0.,0.,0.,1.,1.,1.,1) 
        surface1.intersection(surface2) 
        self.assertTrue(surface.num_vertices()==8 and surface.num_faces()==12 and surface.num_edges()==18) 
        surface2.difference(surface1)        
        self.assertTrue(surface.num_vertices()==14 and surface.num_faces()==24 and surface.num_edges()==36) 


    def test_surface_meshing(self):
        surface = SVMTK.Surface();
        surface.make_sphere(0.0,0.0,0.0,1.0,1.0) 
        self.assertTrue(surface.num_vertices()>0 and surface.num_faces()>0 and surface.num_edges()>0)
        surface.clear()
        surface.implicit_surface(ellipsoid_function, 1.0,10,1.0,1.0 )
        self.assertTrue(surface.num_vertices()>0 and surface.num_faces()>0 and surface.num_edges()>0) 



    def test_span(self):
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,2.,1.,3.,1)       
        self.assertEqual(surface.span(0),(0.,2))
        self.assertEqual(surface.span(1),(0.,1))
        self.assertEqual(surface.span(2),(0.,3))

    def test_fill_holes(self):   
        """ To Be Implemented"""
    def test_triangulate_faces(self):
        """ To Be Implemented"""
    def test_triangulate_hole(self): 
        """ To Be Implemented"""
    def test_surface_clip(self):
        surface =SVMTK.Surface()   
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface.clip(0,0,1.,0,True) 
        self.assertEqual(surface.span(0),(-1.0,1.0))
        self.assertEqual(surface.span(1),(-1.0,1.0))
        self.assertEqual(surface.span(2),(-1.,0.0))



    def test_adjust_boundary(self):
        surface =SVMTK.Surface()   
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface.adjust_boundary(-0.1)
        self.assertTrue(surface.span(0)[1]-surface.span(0)[0] < 2. ) 
        surface.clear() 
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface.adjust_boundary(0.1)
        self.assertTrue(surface.span(0)[1]-surface.span(0)[0] > 2. )


    def test_smoothing(self):
        surface =SVMTK.Surface()   
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface.smooth_laplacian(0.8,1) 
        surface.smooth_taubin(2) 

    def test_mean_curvature_flow(self):
        """ To Be Implemented"""
        
    def test_shortest_surface_path(self):      
        surface =SVMTK.Surface()   
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1)         
        a = surface.shortest_surface_path(SVMTK.Point_3(-1,-1,1), SVMTK.Point_3(1,-1,1))
        #self.assertTrue( a[0]==SVMTK.Point_3(-1,-1,1) and a[1]==SVMTK.Point_3(1,-1,1))

    def test_collapse_and_split_edges(self): 
        surface =SVMTK.Surface()  
        surface.make_cube(0.,0.,0.,1.,1.,1.,1) 
        self.assertTrue(surface.num_vertices()==8 and surface.num_faces()==12 and surface.num_edges()==18) 
        surface.split_edges(0.5)
        self.assertEqual(surface.num_edges(),108) 
        surface.collapse_edges(10)
        self.assertEqual(surface.num_edges(),108) 

    def test_extension(self):
        """ To Be Implemented"""

    def test_reconstruct(self):
        """ To Be Implemented"""

    def test_convex_hull(self):
        """ To Be Implemented"""

    def test_strictly_inside(self):
        surface1=SVMTK.Surface()   
        surface1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface2=SVMTK.Surface()   
        surface2.make_cube(-2.,-2.,-2.,2.,2.,2.,1) 




    def test_seperate_narrow_gaps(self):
        surface=SVMTK.Surface()   
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
    
 


if __name__ == '__main__':
    unittest.main()

