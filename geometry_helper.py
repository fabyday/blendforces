import numpy as np 
import scipy.sparse as sp 

import mesh as mm 


def cotangent_angle_for_vertex(origin_v0, v1, v2):
    # sin / cos = cot :)
    vec1 = v1 - origin_v0
    vec2 = v2 - origin_v0
    return np.cos(vec1, vec2)/np.cross(vec1, vec2)

def make_laplacian(mesh : mm.Mesh):
    # mesh.halfedge.opposite_edge_idx(e_idx)
    vv = mesh.v
    N, dim = vv.shape


    data = []
    rows = [] 
    cols = []
    halfedge = mesh.halfedge 

    def neibour_tri_cot(opposite_e_idx):
        opposite_prev_e_idx = halfedge.prev_edge(opposite_e_idx)
                
        v1_idx, v2_idx = halfedge.edge(opposite_prev_e_idx)
        _, v0_idx = halfedge.edge(opposite_e_idx)
        a_cot = cotangent_angle_for_vertex(vv[v1_idx, :], vv[v2_idx, :], vv[v0_idx, :])
        return a_cot
    
    def cur_tri_cot(e_idx):
        next_edge_idx = halfedge.next_edge(e_idx)
        v0_idx, v1_idx = halfedge.edge(e_idx)
        _, v2_idx = halfedge.edge(next_edge_idx)
        b_cot = cotangent_angle_for_vertex(vv[v1_idx, :], vv[v2_idx, :], vv[v0_idx, :])
        return b_cot 

    for vi in range(N):
        w_sum = 0.0
        e_list = mesh.halfedge.v2e[vi]
        for e_idx in e_list: 
            _, to_vidx = halfedge.edge(e_idx)
            b_cot = cur_tri_cot(e_idx)
            # https://igl.ethz.ch/projects/ARAP/arap_web.pdf

            opposite_e_idx = halfedge.opposite_edge_idx(e_idx)
            a_cot = 0.0
            div_size = 1.0
            if not opposite_e_idx == -1 : 
                a_cot = neibour_tri_cot(opposite_e_idx)
                div_size = 2.0
            w = (b_cot + a_cot)/div_size
            
            rows.append(vi);cols.append(to_vidx);data.append(w)
            w_sum += w

            rows.append(vi); cols.append(vi); data.append(-w_sum)

    
    return sp.csr_matrix((data, (rows, cols)), shape=(N,N), dtype=np.float64)


    