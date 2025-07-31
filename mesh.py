import igl 

import numpy as np 
# import tetgen
import typing 
class Mesh:
    pass 

        

    
class EdgeData:
    def __init__(self,ref : Mesh):
        pass 




class FaceData:
    pass 

class VertexData:
    def __init__(self, v):
        self.__m_v = v 




def append_edge(edge_data : list, edge_meta : dict, edge_idx : int, v1_idx:int, v2_idx : int):
    k = (v1_idx, v2_idx)
    edge_meta[k] = edge_idx
    edge_data[edge_idx] = k

class HalfEdge:
    def __init__(self, v : np.ndarray, f : np.ndarray):
        self.__m_v = v 
        self.__m_f = f



    

    def build(self):
        n_vertices = self.__m_v.shape[0]
        n_faces = self.__m_f.shape[0]
        n_halfedges = n_faces * 3

        # Half-edge arrays
        halfedge_to_vertex = np.zeros(n_halfedges, dtype=int)
        halfedge_to_face = np.zeros(n_halfedges, dtype=int)
        halfedge_next = np.zeros(n_halfedges, dtype=int)
        halfedge_twin = -np.ones(n_halfedges, dtype=int)  # -1 = no twin yet

        face_to_halfedge = np.zeros(n_faces, dtype=int)
        vertex_to_halfedge = -np.ones(n_vertices, dtype=int)  # one outgoing edge per vertex

        # Map from (start, end) to half-edge index (for twin linking)
        edge_map = {}

        for f_id, (i0, i1, i2) in enumerate(self.__m_f):
            # Three half-edge indices for this face
            he0 = f_id * 3 + 0
            he1 = f_id * 3 + 1
            he2 = f_id * 3 + 2

            # Set connections
            halfedge_to_vertex[he0] = i1
            halfedge_to_vertex[he1] = i2
            halfedge_to_vertex[he2] = i0

            halfedge_next[he0] = he1
            halfedge_next[he1] = he2
            halfedge_next[he2] = he0

            halfedge_to_face[he0] = f_id
            halfedge_to_face[he1] = f_id
            halfedge_to_face[he2] = f_id

            face_to_halfedge[f_id] = he0

            # Store outgoing half-edge for vertices
            vertex_to_halfedge[i0] = he0
            vertex_to_halfedge[i1] = he1
            vertex_to_halfedge[i2] = he2

            # For twin linking
            edges = [(i0, i1, he0), (i1, i2, he1), (i2, i0, he2)]
            for start, end, he in edges:
                twin_he = edge_map.get((end, start))
                if twin_he is not None:
                    halfedge_twin[he] = twin_he
                    halfedge_twin[twin_he] = he
                else:
                    edge_map[(start, end)] = he

        
        self.halfedge_to_vertex = halfedge_to_vertex
        self.halfedge_to_face = halfedge_to_face
        self.halfedge_next = halfedge_next
        self.halfedge_twin = halfedge_twin
        self.face_to_halfedge = face_to_halfedge
        self.vertex_to_halfedge = vertex_to_halfedge
    
    def neighbor(self, v_idx):
        he_start =he= self.vertex_to_halfedge[v_idx]
        if he_start == -1:
            return 
        
        faces = []
        neighbors = []

        while True:
            neighbor = self.halfedge_to_vertex[he]
            face = self.halfedge_to_face[he]
            neighbors.append(neighbor)
            faces.append(face)

            # 다음 half-edge로 이동 (twin → next)
            twin = self.halfedge_twin[he]
            if twin == -1:
                break  # boundary reached
            he = self.halfedge_next[twin]
            if he == he_start:
                break  # 순환 종료

        return faces, neighbors 
        
    
    def v2v(self, v_idx):
        return self.__m_v2v[v_idx]

    def e2f(self, e_idx):
        return self.__m_e2f[e_idx]
    
    def prev_edge(self, e_idx):
        return self.__m_e2e[e_idx][0]
    
    def next_edge(self, e_idx):
        return self.__m_e2e[e_idx][-1]

    def v2e(self, v_idx):
        res = [ self.__edge_indexkey_mapper.get([v_idx, to_idx], -1) for to_idx in self.__m_v2v[v_idx] ]
        return list(filter(lambda x : x == -1, res))
    
    def opposite_edge_idx(self, query :typing.Union[ typing.List[int], typing.Tuple[int, int ], int ] ):
        """
        """

        if isinstance(query, tuple) or isinstance(query, list):
            ind = self.__edge_indexkey_mapper.get(query[::-1], -1)

        else : 
            v_indices   = self.__m_edges[query]
            rev_indices = v_indices[::-1]
            ind = self.__edge_indexkey_mapper.get(rev_indices, -1)

        return ind
    

    def edge(self, e_idx):
        """return vertex index"""
        return self.__m_edges[e_idx]
    
    @property
    def edges(self):
        # self.__m_edges
        self.__m_edges = igl.edges(self.__m_f)
        return self.__m_edges



            


class Mesh():
    def __init__(self):
        pass 

    @property
    def v(self):
        return self.__v
    
    @v.setter
    def v(self, vv : np.ndarray):
        self.__v = vv 
    

    def __getitem__(self, idx):
        return self.__v[idx, :]


    @property
    def e(self):
        return self.__m_halfedge.edges
    
    @property
    def f(self):
        return self.__f
    @f.setter
    def f(self, ff : np.ndarray):
        self.__f = ff
    


    def build_halfege(self):
        self.__m_halfedge = HalfEdge(self.__v, self.__f)
        self.__m_halfedge.build()


    
    def load_from_file(self, pth : str):
        v, f = igl.read_triangle_mesh(pth)
        self.__v = v 
        self.__f = f


    @property
    def halfedge(self):
        return self.__m_halfedge

    



# class TetrahedronMesh():
    
#     def load_from_mesh(self, mesh : Mesh):
#         tet = tetgen.TetGen(mesh.v, mesh.f)
#         self.__m_node , self.__m_elem = tet.tetrahedralize(verbose=1)
        

#         print(self.__m_elem)


#     def bulid_individual_aabb_list(self):
#         for elem in self._m_node:
#             elem







if __name__ == "__main__":
    # m = Mesh()
    # m.v = np.random.random((3,3))
    # f =  np.array([[0,1,2]], dtype=np.int32 )
    # m.f =f
    # t=  TetrahedronMesh()
    # t.load_from_mesh(m)


    import igl 

    # v, f = igl.read_triangle_mesh("./data/untitled.obj")
    v, f = igl.read_triangle_mesh("./data/test.obj")
    m = Mesh()
    m.load_from_file("./data/test.obj")
    m.build_halfege()
