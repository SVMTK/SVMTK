
import unittest
import SVMTK


class Domain_Test(unittest.TestCase):

    def test_single_domain(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        domain = SVMTK.Domain(surface_1) 
        self.assertEqual(domain.number_of_surfaces(),1)

    def test_two_domains(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface_2 = SVMTK.Surface() 
        surface_2.make_cube(-2.,-2.,-2.,2.,2.,2.,1) 
        domain = SVMTK.Domain([surface_1,surface_2])
        self.assertEqual(domain.number_of_surfaces(),2)

    def test_mehsing_domains_with_map(self): 
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface_2 = SVMTK.Surface() 
        surface_2.make_cube(-2.,-2.,-2.,2.,2.,2.,1)
        sf= SVMTK.SubdomainMap()
        sf.add("01",3) 
        sf.add("11",2)        
        domain = SVMTK.Domain([surface_1,surface_2],sf)
        domain.create_mesh(1.) 
        self.assertTrue(domain.number_of_cells() >0) 

    def test_get_boundary_and_patches(self): 
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface_2 = SVMTK.Surface() 
        surface_2.make_cube(-2.,-2.,-2.,2.,2.,2.,1)
        sf= SVMTK.SubdomainMap()
        sf.add("01",3) 
        sf.add("11",2)     
        domain = SVMTK.Domain([surface_1,surface_2],sf)
        domain.add_sharp_border_edges(surface_1,0)
        domain.add_sharp_border_edges(surface_2,0)        
        domain.create_mesh(6) 
        surface = domain.get_boundary(2)

        self.assertTrue(surface.num_vertices()==50)
        self.assertTrue(surface.num_faces()==96)
        self.assertTrue(surface.num_edges()==144)  

        surface = domain.get_boundary(1) 

        self.assertEqual(surface.num_vertices(), 0) 
        
        surface = domain.get_boundary(3)
      
        self.assertTrue(surface.num_vertices()==244)
        self.assertTrue(surface.num_faces()==480)
        self.assertTrue(surface.num_edges()==720)  
        
        surfaces =  domain.get_boundaries()  
        self.assertEqual(len(surfaces),2) 

    def test_domain_with_polyline_meshing(self):
        surface = SVMTK.Surface();
        surface.make_sphere(0.0,0.0,0.0,3.0,1.0)
        domain = SVMTK.Domain(surface)
        line0 = [SVMTK.Point_3(0,0,0.0), SVMTK.Point_3(0,1.0,1.0)] 
        domain.add_feature(line0)
        domain.create_mesh(32.)
        self.assertTrue(domain.number_of_curves() >0) 

    def test_domain_meshing(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        domain = SVMTK.Domain(surface_1)
        domain.create_mesh(1)
        self.assertTrue(domain.number_of_cells() >0) 

    def test_domain_meshing_with_borders(self):
        surface = SVMTK.Surface()
        surface.make_cube(0.0,0.0,0.0,1.0,1.0,1.0,1)
        domain = SVMTK.Domain(surface)
        domain.add_sharp_border_edges(surface)
        domain.create_mesh(32.)
        self.assertTrue(domain.number_of_curves() > 0) 

    def test_mesh_lloyd(self):
        surface = SVMTK.Surface()
        surface.make_cube(0.0,0.0,0.0,1.0,1.0,1.0,1.0)
        domain = SVMTK.Domain(surface)
        domain.add_sharp_border_edges(surface)
        domain.create_mesh(6.)
        domain.lloyd() 
       
    def test_mesh_excude(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        domain = SVMTK.Domain(surface_1)
        domain.create_mesh(1)
        domain.exude() 

    def test_mesh_odt(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        domain = SVMTK.Domain(surface_1)
        domain.create_mesh(1)
        domain.odt()        

    def test_mesh_perturb(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        domain = SVMTK.Domain(surface_1)
        domain.create_mesh(1)
        domain.perturb() 
  
    def test_surface_segmentation(self): 
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(0,0,0,1.,1.,1.,1) 
        surface_2 = SVMTK.Surface() 
        surface_2.make_cube(1.,0.,0.,2.,2,2.,1)
        surface_3 = SVMTK.Surface() 
        surface_3.make_cube(0.,1.,0.,2,2.,2.,1)        
        sf= SVMTK.SubdomainMap()
 
        surface_1.union(surface_2)
        surface_1.union(surface_3) 

        domain = SVMTK.Domain(surface_1)
        domain.add_sharp_border_edges(surface_1,85)
        domain.create_mesh(10.) 
        self.assertTrue(domain.number_of_patches()==1)   
        domain.boundary_segmentations()
        self.assertTrue(domain.number_of_patches()==9)       
        



if __name__ == '__main__':
    unittest.main()




