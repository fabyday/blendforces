import primitives as prim

import typing 
import collections 
import numpy as np 
import tetgen
import mesh as mm

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
SPACING_ARG_TYPE = typing.Union[float , str ]
PRIMITIVES_TYPE = typing.Union[mm.Mesh, prim.PrimitiveContainer, typing.List[prim.PrimitiveContainer] ]
MESH_CONVERT_METHOD = str
class SpatialHashing:
    def __init__(self, spacing :  SPACING_ARG_TYPE = "auto", hash_table_size = 100000, default_mesh_converting_method : MESH_CONVERT_METHOD = "AABB" ):
        self.__m_spacing : SPACING_ARG_TYPE = spacing
        self.__m_data  : list = []
        self.__m_hash_table_size = hash_table_size 
        self.__m_hash_table = [ [] for _ in range(hash_table_size) ]
        # self.__m_hierarchy_grid : Grid | None = None

        self.__m_mesh_convert_method = default_mesh_converting_method


        self.__hash_coeff = np.array([73856093, 19349663, 83492791], dtype=np.int32)

    
    def append_primitives(self, primitives:PRIMITIVES_TYPE):
        if isinstance(primitives, list):
            self.__m_data += primitives
        elif isinstance(primitives, prim.PrimitiveContainer) or issubclass(primitives.__class__, prim.PrimitiveContainer) :
            self.__m_data.append(primitives)
        elif isinstance(primitives, mm.Mesh):
            if self.__m_mesh_convert_method == "Tet":
                pass  #TODO
            elif self.__m_mesh_convert_method == "AABB": 
                self.__m_mesh = primitives
                self.__m_data += prim.AABB.mesh_to_primitives(primitives)
    

            
            







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
    
    def __put_primitives_onto_grid_cell(self):
        grid_size = self.__m_grid_size
        grid_idx_list = np.empty((len(self.__m_data), 6))
        for grid_idx, prim in enumerate(self.__m_data) :
            x = (np.floor(prim.__m_minmax/grid_size)).astype(np.int32)
            grid_idx_list [ grid_idx , :  ] = x.reshape(1, -1)





            
            
    def __put_vertice_onto_grid_cell(self):
        grid_size = self.__m_grid_size

        grid_indices = (self.__m_mesh.v / grid_size).astype(np.int32)
        self.__compute_hash_table(grid_indices)

    def __compute_hash(self, xyz):
        """ 
        hash(x,y,z) = ( x p1 xor y p2 xor z p3) mod n
        where p1, p2, p3 are large prime numbers, in
        our case 73856093, 19349663, 83492791, respectively. The value n is the hash table size.
        """
        tmp = np.multiply(xyz, self.__hash_coeff)
        xp1, xp2, xp3 = tmp.ravel()
        return (xp1^xp2^xp3) % self.__m_hash_table_size
    
    def __compute_hash_table(self, idx_list):
        hash_values = [self.__compute_hash(idx_set) for idx_set in idx_list] 
        self.__m_hash_table = {}
        for elem_idx, hash_v in enumerate(hash_values) : 
            value = self.__m_hash_table.get(hash_v, None )
            if value is None :
                value = set()
                self.__m_hash_table[hash_v] = value 
            value.add(elem_idx)


    def query(self, primitives : typing.Union[typing.List[prim.AABB],prim.AABB]):
        if isinstance(primitives, list):
            for prim in primitives:
                self.__intersect_test_phase_1(prim)
                self.__intersect_test_phase_2(prim) 
        else : 
            self.__intersect_test_phase_1(prim)
            self.__intersect_test_phase_2(prim)


    def __intersect_test_phase_1(self, prim : prim.AABB):
        pass #check hashed prims mapped to same hash index that occupied by verts.
        aabb_xyzs = np.floor((prim.data / self.__m_grid_size)).astype(np.int32)
        aabb_hash = [self.__compute_hash(aabb) for aabb in aabb_xyzs]
        hashes = self.__m_hash_table[aabb_hash]



        return 




    def __intersect_test_phase_2(self, prim, v_idx):
        pass 

    def precompute(self):
        self.__m_grid_size = self.__compute_grid_size()
        self.__put_vertice_onto_grid_cell()





    
if __name__ == "__main__":

    import os , glob 
    data_path = "D:\\lab\\2022\\mycode\\FaceCaptureWithIK\\data\\ICT-data"
    neutral_pth = os.path.join(data_path, "generic_neutral_mesh.obj")
    shapes_path = os.path.join(data_path, "shapes")
    file_pths = glob.glob(os.path.join(shapes_path, "**.obj"))
    neutral = mm.Mesh()
    neutral.load_from_file(neutral_pth)
    sh = SpatialHashing()
    sh.append_primitives(neutral)
    sh.precompute()


    sh.query()