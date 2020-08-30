import SVMTK as svm


if __name__ == "__main__":

   surf1 = svm.Surface()

   surf1.make_cube(0,0,0,2,2,2) 

   maker = svm.Domain(surf1)

   maker.add_sharp_border_edges(surf1,70)

   maker.create_mesh(126) # need high resolution

   maker.save("box.mesh")

