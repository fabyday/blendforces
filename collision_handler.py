import primitives as prim

import typing 
import collections 
import numpy as np 
import mesh as mm

# class Grid():
#     def __init__(self, grid_length):
#         self.__m_grid_length = grid_length 
#         self.__m_grid = None
#         self

#     def add(self, primitives):
#         if isinstance(primitives, list):
#             self.__m_data_list += primitives
#         elif isinstance(primitives, prim.PrimitiveContainer) or issubclass(primitives, prim.PrimitiveContainer) :
#             self.__m_data_list.append(primitives)

#     def build(self):
#         for prim_data  in self.__m_data:
#             np.floor(prim_data.data / self.__m_grid_size)





class CollisionEnergy:
    """
        impl https://la.disneyresearch.com/wp-content/uploads/Robust-Treatment-of-Simultaneous-Collisions-Paper.pdf
    """
    @staticmethod
    def momentum(Mass, velocity_t : np.ndarray):
        """
            M : diagonal mass matrix.
            velocity : velocity of vertex positions. shape[N,3]

        """
        p_t = Mass @ velocity_t.reshape(-1,3)

        return p_t.reshape(-1,3)
    
    @staticmethod
    def kinetic_energy(mass, momentum_engery_t):
        """
            M : spase diagonal mass matrix. [N, N]
            momentum engery at t time. N,3
        """
        flat_momentum = momentum_engery_t.reshape(1,-1)
        mass_dot_momentum = ( (1.0 / mass) @ momentum_engery_t).reshape(-1,1)
        return flat_momentum@mass_dot_momentum
    
    @staticmethod
    def tri_vert_collision_energy(face_v_idx, face_normal, vertex_ind, vv):
        """
            if return 1, then it is not contatct 
            if return 0, then it on the exact triangle plane.
            if return -1, then it is penetrating state.
            see also http://wscg.zcu.cz/wscg1996/papers96/Zacharias_96.pdf
        """
        v1v2v3 =  vv[face_v_idx, : ]
        v1 = v1v2v3[0, :]
        v2 = v1v2v3[1, :]
        v3 = v1v2v3[2, :]
        # plane equation n@x - bias(:=n@v) = 0 
        bias = v1@face_normal
        v4 = vv[vertex_ind]

        # project on triangle surface 
        # see also https://math.stackexchange.com/questions/100761/how-do-i-find-the-projection-of-a-point-onto-a-plane
        target_v = (v4 - v1)
        projected_v3 = v1 + target_v - target_v @ face_normal / np.linalg.normI(face_normal)


        def tri_size(e1, e2):
            e2_norm = np.linalg.norm(e2)
            h = e1 - (e1.T @ e2)@e2/(e2_norm**2)
            tri_area = h*e2_norm*0.5
            return tri_area            
        

        tri_area = tri_size(v2-v1,v3-v1)
        bcp = tri_size(v2-v4,v3-v4)/tri_area # lambda 1
        acp = tri_size(v3-v4,v1-v4)/tri_area # lambda 2 
        abp = tri_size(v1-v4,v2-v4)/tri_area # lambda 3


        if bcp + acp + abp > 1:
            "be-point was out of triangle"
        







        
        barecentric =  np.array([bcp, acp, abp]).reshape(-1,1)

        return face_normal @ (v4 - v1v2v3@barecentric)
        

    @staticmethod
    def edge_edge_collision_energy(edge1_indice, edge2_indice, vv):
        # see also https://paulbourke.net/geometry/pointlineplane/
        e11i, e12i = edge1_indice
        e21i, e22i = edge2_indice
        v11 = vv[e11i]
        v12 = vv[e12i]
        v21 = vv[e21i]
        v22 = vv[e22i]

        e1_dir = v12 - v11
        e1_norm_dir = e1_dir / np.linalg.norm(e1_dir)
        e2_dir = v22 - v21
        e2_norm_dir = e2_dir/np.linalg.norm(e2_dir)
        
        
        
        (v11-v21)


        

    


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


    def __intersect_test_phase_1(self, primitive : prim.AABB):
        pass #check hashed prims mapped to same hash index that occupied by verts.
        aabb_xyzs = np.floor((primitive.data / self.__m_grid_size)).astype(np.int32)
        min_ = aabb_xyzs[prim.AABB.MIN_IDX]
        max_ = aabb_xyzs[prim.AABB.MAX_IDX]
        hash_list = []
        for i in range(min_[..., 0], max_[..., 0] + 1)        :
            for j in range(min_[..., 1], max_[..., 1] + 1)        :
                for k in range(min_[..., 2], max_[..., 2] + 1)        : 
                    hh = np.array([i,j,k])
                    hash_list.append(hh)

        # aabb_hash = [self.__compute_hash(aabb) for aabb in aabb_xyzs]
        aabb_hash = [self.__compute_hash(aabb) for aabb in hash_list]
        col_candidate_v_indice = []
        for cel_hash in aabb_hash:
            idx = self.__m_hash_table[aabb_hash]
            col_candidate_v_indice.append(col_candidate_v_indice)
        return col_candidate_v_indice
        


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




    