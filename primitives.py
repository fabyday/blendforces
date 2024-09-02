
import numpy as np 
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
    def __init__(self, copy=False ):
        super().__init__(copy=copy)
        self.__m_minmax = np.empty(2,3, dtype=np.float64)
    


    def __build(self):
        pass

    def max_length(self):
        length = np.abs(np.subtract(self.__m_minmax[0, :], self.__m_minmax[1, :]))
        return np.max(length)

        


class Sphere(PrimitiveContainer):
    def __init__(self, copy=False):
        super().__init__(copy=copy)
        pass 