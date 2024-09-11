
import numpy as np 
import mesh as mm
import tetgen
class PrimitiveContainer:
    def __init__(self, copy = False ):
        self.__m_data = None 
        self.__m_copy = copy 
    
    

    def __build(self):
        return NotImplemented


    @property
    def data(self):
        return self.__m_data
    
    @data.setter
    def data(self, data):
        *_, dim = data.shape
        
        if not dim == 3 : 
            raise ValueError("dimension is not 3-D.")
        
        if self.__m_copy : 
            self.__m_data = np.copy(data)
        else :
            self.__m_data = data 
        

    

    @staticmethod
    def mesh_to_primitives(mesh):
        return 



class Tetrahedron(PrimitiveContainer):
    def __init__(self, copy = False ):
        super().__init__(copy=copy)


    def __build(self):
        n, dim = self.__m_data.shape
        if n  == 3 : 
            tet = tetgen.TetGen(self.__m_data, [0,1,2])
            self.__m_node, self.__m_elem = tet.tetrahedralize()
        else :
            raise ValueError("shape was not 3 vertex.")




class AABB(PrimitiveContainer):
    MIN_IDX = 0
    MAX_IDX = 1
    def __init__(self, copy=False ):
        super().__init__(copy=copy)
        self.__m_minmax = np.empty((2,3), dtype=np.float64)
        self.__m_meta = {}


    def __build(self):
        self.__m_minmax[AABB.MIN_IDX, :] = np.min(self.data, axis=0)
        self.__m_minmax[AABB.MAX_IDX, :] = np.max(self.data, axis=0)
        

    def max_length(self):
        length = np.abs(np.subtract(self.__m_minmax[0, :], self.__m_minmax[1, :]))
        return np.max(length)


    def set(self, key, value):
        self.__m_meta[key] = value
    def get(self, key):
        return self.__m_meta[key]

    @staticmethod
    def mesh_to_primitives(mesh : mm.Mesh , span_length= 0.01):
        vvf = mesh.v[mesh.f.reshape(-1), : ]
        ffv = vvf.reshape(-1, 9)
        prims = []
        for v3, ff in zip(ffv, mesh.f):
            vv3 = v3.reshape(3,3)
            aabb = AABB()
            aabb.data = vv3 
            aabb.__build()
            aabb.set("f", ff.reshape(-1)) # add meta
            prims.append(aabb)
        return prims


        


class Sphere(PrimitiveContainer):
    def __init__(self, copy=False):
        super().__init__(copy=copy)
        pass 