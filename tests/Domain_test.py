
import unittest
import SVMTK
import os

tests_dir = os.path.dirname(__file__)

class Domain_Test(unittest.TestCase):

    def test_unpack_mesh(self):
        domain = SVMTK.Domain(f"{tests_dir}/Data/cube.mesh") 
        
        points = domain.get_points()
        cells  = domain.get_cells()
        facets = domain.get_facets(False)
        marked = domain.get_facets(True)

        self.assertEqual(points.shape, (domain.num_vertices(),3))        
        # CGAL doesn't load all facets and cells
        del domain

    def test_refine(self):
        domain = SVMTK.Domain(f"{tests_dir}/Data/cube.mesh") 
        num_cells = domain.num_cells()
        domain.refine(0.5)
        self.assertTrue(domain.num_cells()>num_cells)
        del domain 

    def test_single_domain(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        domain = SVMTK.Domain(surface_1) 
        self.assertEqual(domain.num_surfaces(),1)
        del domain

    def test_two_domains(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        surface_2 = SVMTK.Surface() 
        surface_2.make_cube(-2.,-2.,-2.,2.,2.,2.,1) 
        domain = SVMTK.Domain([surface_1,surface_2])
        self.assertEqual(domain.num_surfaces(),2)
        del domain

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
        self.assertTrue(domain.num_cells() >0) 
        del domain

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
        domain.create_mesh(1) 
        surface = domain.get_boundary(2)
        
        self.assertTrue(surface.num_vertices()>0)
        self.assertTrue(surface.num_faces()>0)
        self.assertTrue(surface.num_edges()>0)  
        
        surface = domain.get_boundary(1) 
        
        self.assertEqual(surface.num_vertices(), 0) 
        
        surface = domain.get_boundary(3)
        
        self.assertTrue(surface.num_vertices()>0)
        self.assertTrue(surface.num_faces()>0)
        self.assertTrue(surface.num_edges()>0)  
        
        surfaces =  domain.get_boundaries()  
        self.assertEqual(len(surfaces),2) 
        del domain

    def test_domain_with_polyline_meshing(self):
        surface = SVMTK.Surface();
        surface.make_sphere(0.0,0.0,0.0,3.0,1.0)
        domain = SVMTK.Domain(surface)
        line0 = [SVMTK.Point_3(0,0,0.0), SVMTK.Point_3(0,1.0,1.0)] 
        domain.add_feature(line0)
        domain.create_mesh(1.0)
        self.assertTrue(domain.num_curves() >0) 
        del domain

    def test_domain_meshing(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        domain = SVMTK.Domain(surface_1)
        domain.create_mesh(1)
        self.assertTrue(domain.num_cells() >0) 
        del domain


    def test_domain_meshing_with_borders(self):
        surface = SVMTK.Surface()
        surface.make_cube(0.0,0.0,0.0,1.0,1.0,1.0,1.0)
        domain = SVMTK.Domain(surface)
        domain.add_sharp_border_edges(surface,60)
        domain.create_mesh(1)
        self.assertTrue(domain.num_curves() > 0) 
        del domain

    def test_mesh_lloyd(self):
        surface = SVMTK.Surface()
        surface.make_sphere(0.0,0.0,0.0,3.0,1.0)
        domain = SVMTK.Domain(surface)
        domain.create_mesh(1.0)
        domain.lloyd(time_limit=1.0) 
        del domain 
       
       
    def test_mesh_excude(self):
        surface = SVMTK.Surface() 
        surface.make_cube(-1.,-1.,-1.,1.,1.,1.,1.0) 
        domain = SVMTK.Domain(surface)
        domain.add_sharp_border_edges(surface)
        domain.create_mesh(1.0) 
        domain.exude() 
        del domain

    def test_mesh_odt(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        domain = SVMTK.Domain(surface_1)
        domain.create_mesh(1)
        domain.odt()        
        del domain 

    def test_mesh_perturb(self):
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(-1.,-1.,-1.,1.,1.,1.,1) 
        domain = SVMTK.Domain(surface_1)
        domain.create_mesh(1)
        domain.perturb() 
        del domain
  
    def test_surface_segmentation(self): 
        surface_1 = SVMTK.Surface() 
        surface_1.make_cube(0,0,0,1.,1.,1.,1.0) 
        surface_2 = SVMTK.Surface() 
        surface_2.make_cube(1.,0.,0.,2.,2,2.,1.0)
        surface_3 = SVMTK.Surface() 
        surface_3.make_cube(0.,1.,0.,2,2.,2.,1.0)        
        sf= SVMTK.SubdomainMap()
 
        surface_1.union(surface_2)
        surface_1.union(surface_3) 

        domain = SVMTK.Domain(surface_1)
        domain.add_sharp_border_edges(surface_1,85)
        domain.create_mesh(1.0) 
        self.assertTrue(domain.num_patches()==1)   
        domain.boundary_segmentations()
        self.assertTrue(domain.num_patches()==9)     
          
        del domain 



if __name__ == '__main__':
    unittest.main()




