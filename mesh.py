import igl 

import numpy as np 
import tetgen
class Mesh():
    def __init__(self):
        pass 



    @property
    def v(self):
        return self.__v
    
    @v.setter
    def v(self, vv : np.ndarray):
        self.__v = vv 
    @property
    def f(self):
        return self.__f
    @f.setter
    def f(self, ff : np.ndarray):
        self.__f = ff
 

    def load_from_file(self, pth : str):
        v, f = igl.read_triangle_meshes(pth)
        self.__v = v 
        self.__f = f


class TetrahedronMesh():
    
    def load_from_mesh(self, mesh : Mesh):
        tet = tetgen.TetGen(mesh.v, mesh.f)
        self.__m_node , self.__m_elem = tet.tetrahedralize(verbose=1)
        

        print(self.__m_elem)


    def bulid_individual_aabb_list(self):
        for elem in self._m_node:
            elem







if __name__ == "__main__":
    # m = Mesh()
    # m.v = np.random.random((3,3))
    # f =  np.array([[0,1,2]], dtype=np.int32 )
    # m.f =f
    # t=  TetrahedronMesh()
    # t.load_from_mesh(m)
    v = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
     
    ]
    )
    f = np.vstack(
        [
            [0, 1, 2],
        
        ]
    )
    tgen = tetgen.TetGen(v, f)
    nodes, elems = tgen.tetrahedralize()
    tgen.grid.plot(show_edges=True)