from . import primitives as prim

import typing 
import collections 
import numpy as np 
import tetgen

class Grid():
    def __init__(self, grid_length):
        self.__m_grid_length = grid_length 
        self.__m_grid = None
        self

    def add(self, primitives):
        if isinstance(primitives, list):
            self.__m_data_list += primitives
        elif isinstance(primitives, prim.PrimitiveContainer) or issubclass(primitives, prim.PrimitiveContainer) :
            self.__m_data_list.append(primitives)

    def build(self):
        for prim_data  in self.__m_data:
            np.floor(prim_data.data / self.__m_grid_size)


        
    



# impl Hierarchical spatialhashing according to paper https://matthias-research.github.io/pages/publications/tetraederCollision.pdf.
SPACING_ARG_TYPE = typing.Union[float , typing.Literal["auto"] ]

class SpatialHashing:
    def __init__(self, spacing :  SPACING_ARG_TYPE = "auto", hash_table_size = 10000):
        self.__m_spacing : SPACING_ARG_TYPE = spacing
        self.__m_data  : list = []
        self.__m_hash_table_size = hash_table_size 
        self.__m_hash_table_size = [ [] for _ in hash_table_size ]
        self.__m_hierarchy_grid : Grid | None = None


        self.__hash_coeff = np.array([73856093, 19349663, 83492791], dtype=np.uint32)

    
    def append_primitives(self, primitives :prim.PrimitiveContainer|list[prim.PrimitiveContainer]):
        if isinstance(primitives, list):
            self.__m_data += primitives
        elif isinstance(primitives, prim.PrimitiveContainer) or issubclass(primitives, prim.PrimitiveContainer) :
            self.__m_data.append(primitives)






    def __compute_grid_size(self):
        if self.__m_spacing == "auto":

            num_data = len(self.__m_data)
            avg_l = 0.0
            for prim in self.__m_data:
                avg_l += prim.max_length()
            avg_l /= num_data 

            grid_size = avg_l

        else : 
            grid_size = self.__m_spacing
            
        return grid_size
    
    def __put_primitives_onto_grid_cell(self, grid_size):
        pass 



    def __compute_hash(self, xyz, hash_table_size):
        """ 
        hash(x,y,z) = ( x p1 xor y p2 xor z p3) mod n
        where p1, p2, p3 are large prime numbers, in
        our case 73856093, 19349663, 83492791, respectively. The value n is the hash table size.
        """
        tmp = np.multiply(xyz, self.__hash_coeff)
        xp1, xp2, xp3 = tmp.ravel()
        return (xp1^xp2^xp3) % hash_table_size
    
    def __compute_hash_table(self, hash_table_size):
        
        self.__m_data

    def precompute(self):

        self.__m_grid_size = self.__compute_grid_size()
        self.__put_primitives_onto_grid_cell(self.__m_grid_size)
        self.__compute_hash_table(self.__m_hash_table_size)














    



    