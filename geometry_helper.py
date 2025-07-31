import numpy as np 
import scipy.sparse as sp 

import mesh as mm 


import igl 
def compute_vertex_normals(vertices, faces):
    vertex_normals = np.zeros_like(vertices)

    for face in faces:
        i0, i1, i2 = face
        v0, v1, v2 = vertices[i0], vertices[i1], vertices[i2]

        # 면의 normal 계산 (삼각형의 두 변의 외적)
        edge1 = v1 - v0
        edge2 = v2 - v0
        face_normal = np.cross(edge1, edge2)
        face_normal = face_normal / (np.linalg.norm(face_normal) + 1e-10)  # 정규화

        # 정점에 더해줌
        vertex_normals[i0] += face_normal
        vertex_normals[i1] += face_normal
        vertex_normals[i2] += face_normal

    # 정점별로 정규화
    vertex_normals = np.array([
        n / (np.linalg.norm(n) + 1e-10)
        for n in vertex_normals
    ])

    return vertex_normals



def cotangent_angle_for_vertex(origin_v0, v1, v2, eps=1e-8):
    # sin / cos = cot :)
    vec1 = v1 - origin_v0
    vec2 = v2 - origin_v0

    norm1 = np.linalg.norm(vec1)
    norm2 = np.linalg.norm(vec2)

    if norm1 < eps or norm2 < eps:
        return 0.0  # 거의 겹친 점 처리

    vec1 /= norm1
    vec2 /= norm2

    cos = np.dot(vec1, vec2)
    sin = np.linalg.norm(np.cross(vec1, vec2))

    if sin < eps:
        return 0.0  # 각도가 너무 작아서 수치적으로 위험

    return cos / sin

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

    
    return sp.csc_matrix((data, (rows, cols)), shape=(N,N), dtype=np.float64)

def make_neibor_matrix(mesh : mm.Mesh):
    
    rows = []
    cols = [] 
    datas = []
    edge_nums = 0
    neutral = mesh
    w_list = [] 
    for row, v_i in enumerate(range(len(mesh.v))):
        f, edge_list = mesh.halfedge.neighbor(v_i)
        for ei in range(len(edge_list)):
            if(ei == - 1):
                continue 
            prev_e_idx = ei - 1
            next_e_idx = (ei + 1) % len(edge_list)
            prev_v = edge_list[prev_e_idx]
            cur_v = edge_list[ei]
            next_v = edge_list[next_e_idx]
            wij = 0 
            
            if(prev_v != -1):
                wij+=cotangent_angle_for_vertex(neutral.v[prev_v], neutral.v[v_i], neutral.v[cur_v] )
            if(next_v != - 1):
                wij+=cotangent_angle_for_vertex(neutral.v[next_v], neutral.v[v_i], neutral.v[cur_v] )
            wij *= 0.5
            w_list.append(wij)
           
        for v_j in edge_list:
            edge_nums += 1
            rows.append(row)
            rows.append(row)
            cols.append(v_i)
            cols.append(v_j)
            datas.append(1)
            datas.append(-1)
            
    return np.array(w_list), sp.coo_matrix((datas, (rows, cols)), shape=( edge_nums , len(mesh.v)) ).tocsc()

def make_laplacian(mesh : mm.Mesh):
    """
        COTANGENT MAT
    """
    v= mesh.v 
    f = mesh.f
    # mesh.halfedge.v2e(v)
    L = -igl.cotmatrix(v,f)
    eps = 1e-6 # reg term
    L += eps * sp.eye(L.shape[0])  # εI 더해서 definite하게 만듬

    return L
    print(L[0,0])
    print(L[1,1])
    print(L[-1,-1])
    G = igl.grad(v,f)
    dblA = igl.doublearea(v,f)
    
    M = igl.massmatrix(v, f, igl.MASSMATRIX_TYPE_VORONOI)
    from scipy.sparse.linalg import inv

    # 정규화된 Cotangent Laplacian (optional)
    L_norm = inv(M) @ L
    return L_norm
    # return igl.cotmatrix(v,f)
import matplotlib.pyplot as plt
import scipy.sparse as sp
import numpy as np

def visualize_signed_sparse_matrix(A, title="Signed Sparse Matrix"):
    if not sp.isspmatrix(A):
        raise ValueError("Input must be a scipy sparse matrix.")

    A = A.tocoo()  # COO 형식으로 변환

    # 색상 매핑: 양수는 빨강, 음수는 파랑
    colors = np.where(A.data > 0, 'red', 'blue')

    plt.figure(figsize=(6, 6))
    plt.scatter(A.col, A.row, c=colors, s=5)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("Columns")
    plt.ylabel("Rows")
    plt.grid(False)
    plt.show()
if __name__ == "__main__":
    import os , glob 
    data_path = "D:\\lab\\2022\\mycode\\FaceCaptureWithIK\\data\\ICT-data"
    neutralpth = os.path.join(data_path, "generic_neutral_mesh.obj")
    neutral = mm.Mesh()
    neutral.load_from_file(neutralpth)
    print(type(make_laplacian(neutral)))






