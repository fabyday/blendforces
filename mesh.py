import igl 

import numpy as np 
import tetgen
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
        self.__m_edges = [[] for _ in range(3* len(self.__m_f))] 
        self.__edge_indexkey_mapper = {}
        self.__m_e2f = [[] for _ in range(len(self.__m_edges))] 
        self.__m_e2e = [[] for _ in range(len(self.__m_edges))]

        self.__m_v2f =  [ [] for _ in range(len(self.__m_v))]
        self.__m_v2v = [set() for _ in range(len(self.__m_v))]


        self.f2v = self.__m_f


        # self.__m_edges = igl.edges(self.__m_f)


        for fi, (i, j, k) in enumerate(self.__m_f):
            e_offset = fi
            
            append_edge(self.__m_edges, self.__edge_indexkey_mapper, e_offset*3, i, j)
            append_edge(self.__m_edges, self.__edge_indexkey_mapper, e_offset*3 + 1, j, k)
            append_edge(self.__m_edges, self.__edge_indexkey_mapper, e_offset*3 + 2, k, i)
            
            self.__m_e2e[e_offset*3] += [-1, e_offset*3+1]
            self.__m_e2e[e_offset*3+1] += [e_offset*3 , e_offset*3 +2]
            self.__m_e2e[e_offset*3+2] += [e_offset *3+1 , -1]
            
            self.__m_e2f[ e_offset*3   ] = fi
            self.__m_e2f[ e_offset*3+1 ] = fi
            self.__m_e2f[ e_offset*3+2 ] = fi



            self.__m_v2v[i] += [j,k]
            self.__m_v2v[j] += [k,i]
            self.__m_v2v[k] += [i,j]


            self.__m_v2f[i].append(fi)  
            self.__m_v2f[j].append(fi)
            self.__m_v2f[k].append(fi)

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
    def e(self):
        pass
    
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


    import igl 

    # v, f = igl.read_triangle_mesh("./data/untitled.obj")
    v, f = igl.read_triangle_mesh("./data/test.obj")
    m = Mesh()
    m.load_from_file("./data/test.obj")
    m.build_halfege()
