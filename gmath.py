from numba import njit 
import numpy as np 
import scipy.sparse as sp 
import geometry_helper as geo 
import mesh 
def make_sparse_matrix_triplet_function(raw_list, col_list, data_list):
    def append_data(i, j, val):
        raw_list.append(3*i + 0); col_list.append(3*j + 0) ; data_list.append(val)
        raw_list.append(3*i + 1); col_list.append(3*j + 1) ; data_list.append(val)
        raw_list.append(3*i + 2); col_list.append(3*j + 2) ; data_list.append(val)
    return append_data


def bend_constraint_A(neutral : mesh.Mesh, kb : float):
    vsize = rowsize = len(neutral.v)
    dim = 3
    rows = []
    cols = [] 
    datas = [] 
    def append_block(x,row_offset, col_offset):
        nonlocal rows, cols, datas
        rows.append(row_offset+0);rows.append(row_offset+1);rows.append(row_offset+2)
            
        cols.append(col_offset+0);cols.append(col_offset+1);cols.append(col_offset+2)
        datas.append(x);datas.append(x);datas.append(x)
        
    # xi, row_num is same, For convenience, I aliased it.
    for row_num, xi in enumerate(range(vsize)):
        _, edges = neutral.halfedge.neighbor(xi)
        for ei in range(len(edges)):
            prev_e_idx = ei - 1
            next_e_idx = (ei + 1) % len(edges)
            prev_vj = edges[prev_e_idx]
            cur_vj = edges[ei]
            next_vj = edges[next_e_idx]
            wij = 0
            if(prev_vj != -1):
                wij+=geo.cotangent_angle_for_vertex(neutral.v[prev_vj], neutral.v[xi], neutral.v[cur_vj] )
            if(next_vj != - 1):
                wij+=geo.cotangent_angle_for_vertex(neutral.v[next_vj], neutral.v[xi], neutral.v[cur_vj] )
            wij *= 0.5
            # e = neutral.v[xi] - neutral.v[cur_vj]
            
            block_row_offset = row_num * dim
            block_col1_offset  = xi * dim  # xi*3
            block_col2_offset  = cur_vj * dim  # j*3
            # block shape
            append_block(wij, block_row_offset, block_col1_offset)
            append_block(-wij, block_row_offset, block_col2_offset)
            
        return sp.coo_matrix((datas, (rows, cols)), shape=(3*vsize, 3*vsize)).tocsc()
    # data = []
    # row = []
    # col = []

    # L_coo = L.tocoo()

    # for i, j, v in zip(L_coo.row, L_coo.col, L_coo.data):
    #     for axis in range(3):
    #         row.append(3 * i + axis)
    #         col.append(3 * j + axis)
    #         data.append(v)

    # return sp.coo_matrix((data, (row, col)), shape=(R, C)).tocsc()
    

@njit
def _bend_constraint_b_dense_part(R_list, dvv):
    """differential coordinates : dvv"""
    dvvv = dvv.reshape(-1,3)
    res = np.zeros_like(dvvv)
    
    for i in range(R_list.shape[0]):
        a = R_list[i] @ dvvv[i, :].reshape(-1,1)
        res[i, :] = a.T
    return res 

@njit
def _bend_constraint_b_get_Rt(S):
    reshaped_S = S.reshape(-1, 3, 3)
    res = np.zeros_like(reshaped_S)
    for i in range(reshaped_S.shape[0]):
        Si = reshaped_S[i, ...]
        u, s, vh = np.linalg.svd(Si)
        R = vh.T @ u.T
        if np.linalg.det(R) < 0:
            vh[-1, :] *= -1
            # u[-1, :] *= -1
            # R = vh @ u.T
            R = vh.T @ u.T
        res[i,:, :] = R
        

    
    return res

def linearize_arap_params(neutral : mesh.Mesh):
    vsize = rowsize = len(neutral.v)
    block_row = 9
    block_col = 3
    dim = 3
    rows = []
    cols = [] 
    datas = [] 
    def append_block(x,row_offset, col_offset):
        nonlocal rows, cols, datas
        rows.append(row_offset+0);rows.append(row_offset+1);rows.append(row_offset+2)
        rows.append(row_offset+3);rows.append(row_offset+4);rows.append(row_offset+5)
        rows.append(row_offset+6);rows.append(row_offset+7);rows.append(row_offset+8)
            
        cols.append(col_offset+0);cols.append(col_offset+1);cols.append(col_offset+2)
        cols.append(col_offset+0);cols.append(col_offset+1);cols.append(col_offset+2)
        cols.append(col_offset+0);cols.append(col_offset+1);cols.append(col_offset+2)
        datas.append(x[0]);datas.append(x[0]);datas.append(x[0])
        datas.append(x[1]);datas.append(x[1]);datas.append(x[1])
        datas.append(x[2]);datas.append(x[2]);datas.append(x[2])
        
    # xi, row_num is same, For convenience, I aliased it.
    for row_num, xi in enumerate(range(vsize)):
        _, edges = neutral.halfedge.neighbor(xi)
        for ei in range(len(edges)):
            prev_e_idx = ei - 1
            next_e_idx = (ei + 1) % len(edges)
            prev_vj = edges[prev_e_idx]
            cur_vj = edges[ei]
            next_vj = edges[next_e_idx]
            wij = 0
            if(prev_vj != -1):
                wij+=geo.cotangent_angle_for_vertex(neutral.v[prev_vj], neutral.v[xi], neutral.v[cur_vj] )
            if(next_vj != - 1):
                wij+=geo.cotangent_angle_for_vertex(neutral.v[next_vj], neutral.v[xi], neutral.v[cur_vj] )
            wij *= 0.5
            e = neutral.v[xi] - neutral.v[cur_vj]
            
            block_row_offset = row_num * 9
            block_col1_offset  = xi * dim  # xi*3
            block_col2_offset  = cur_vj * dim  # j*3
            # block shape
            #[[x  ]......[other blocks...]]
            #[[ x ]......[other blocks...]]
            #[[  x]......[other blocks...]]
            #[[y  ]......[other blocks...]]
            #[[ y ]......[other blocks...]]
            #[[  y]......[other blocks...]]
            #[[z  ]......[other blocks...]]
            #[[ z ]......[other blocks...]]
            #[[  z]......[other blocks...]]
            append_block(wij*e, block_row_offset, block_col1_offset)
            append_block(-wij*e, block_row_offset, block_col2_offset)
            
        return sp.coo_matrix((datas, (rows, cols)), shape=(9*vsize, 3*vsize)).tocsc()
            
def precompute_arap_synbolic_S_sums( neutral : mesh.Mesh):
    vsize = rowsize = len(neutral.v)
    datas = [] 
    rows = []
    cols = []
    dim = 3
    R_mat_size = 9
    row_size = dim*vsize
    col_size = R_mat_size * vsize # 9 is R matrix size. 3x3
    def append_block(x,row_offset, col_offset):
        nonlocal rows, cols, datas
        
        rows.append(row_offset+0);rows.append(row_offset+0);rows.append(row_offset+0)
        rows.append(row_offset+1);rows.append(row_offset+1);rows.append(row_offset+1)
        rows.append(row_offset+2);rows.append(row_offset+2);rows.append(row_offset+2)
            
        cols.append(col_offset+0);cols.append(col_offset+1);cols.append(col_offset+2)
        cols.append(col_offset+3);cols.append(col_offset+4);cols.append(col_offset+5)
        cols.append(col_offset+6);cols.append(col_offset+7);cols.append(col_offset+8)
        datas.append(x[0]);datas.append(x[1]);datas.append(x[2])
        datas.append(x[0]);datas.append(x[1]);datas.append(x[2])
        datas.append(x[0]);datas.append(x[1]);datas.append(x[2])
        # datas.append(x[0]);datas.append(x[0]);datas.append(x[0])
        # datas.append(x[1]);datas.append(x[1]);datas.append(x[1])
        # datas.append(x[2]);datas.append(x[2]);datas.append(x[2])
    for row_num, xi in enumerate(range(vsize)):
        _, edges = neutral.halfedge.neighbor(xi)
        for ei in range(len(edges)):
            prev_e_idx = ei - 1
            next_e_idx = (ei + 1) % len(edges)
            prev_vj = edges[prev_e_idx]
            cur_vj = edges[ei]
            next_vj = edges[next_e_idx]
            wij = 0
            if(prev_vj != -1):
                wij+=geo.cotangent_angle_for_vertex(neutral.v[prev_vj], neutral.v[xi], neutral.v[cur_vj] )
            if(next_vj != - 1):
                wij+=geo.cotangent_angle_for_vertex(neutral.v[next_vj], neutral.v[xi], neutral.v[cur_vj] )
            wij *= 0.5
            # e = neutral.v[cur_vj] - neutral.v[xi]
            e = neutral.v[xi]-neutral.v[cur_vj]

            block_row_offset = row_num * dim
            block_col1_offset  = xi * R_mat_size  # xi*9
            block_col2_offset  = cur_vj * R_mat_size  # j*9
            value = wij*e*0.5
            append_block(value, block_row_offset, block_col1_offset)
            append_block(value, block_row_offset, block_col2_offset)
            
        return sp.coo_matrix((datas, (rows, cols)), shape=(row_size, col_size)).tocsc()
    
def bend_constraint_b(b : np.ndarray, xt : np.ndarray, L : sp.csc_matrix, neutral : mesh.Mesh,  precomputed_rhs : sp.csc_matrix, precomtued_neighbor_mat : sp.csc_matrix):
# def bend_constraint_b(b : np.ndarray, xt : np.ndarray, L : sp.csc_matrix, neutral : mesh.Mesh, kb : float, precomputed_rhs : sp.csc_matrix, precomtued_neighbor_mat : sp.csc_matrix):
    n_pose = neutral.v.reshape(1,-1)
    col = []
    row = []
    data = []
    S = (precomtued_neighbor_mat@xt.reshape(-1,1)).reshape(-1,3,3)
    
    # # S2 = precomuted_w_edges[:,:,np.newaxis] @ (precomtued_neighbor_mat@xt)[:, np.newaxis, :]
    # S = np.zeros((len(xt), 3,3))
    # for i in range(len(xt)):
    #     faces, edges = neutral.halfedge.neighbor(i)
    
    #     # ssss = neutral.f[faces]
    #     # covmat = np.zeros((3,3))
    #     for ei in range(len(edges)):
    #         if(ei == - 1):
    #             continue 
    #         prev_e_idx = ei - 1
    #         next_e_idx = (ei + 1) % len(edges)
    #         prev_v = edges[prev_e_idx]
    #         cur_v = edges[ei]
    #         next_v = edges[next_e_idx]
    #         wij = 0 
            
    #         if(prev_v != -1):
    #             wij+=geo.cotangent_angle_for_vertex(neutral.v[prev_v], neutral.v[i], neutral.v[cur_v] )
    #         if(next_v != - 1):
    #             wij+=geo.cotangent_angle_for_vertex(neutral.v[next_v], neutral.v[i], neutral.v[cur_v] )
    #         wij *= 0.5
    #         e = neutral.v[cur_v] - neutral.v[i]
    #         e_p = xt[cur_v] - xt[i]
    #         S[i, :, :] += wij*e.reshape(-1,1) @ e_p.reshape(1,-1)
            
    # Lxt = L @ xt.reshape(-1,1)
    # S = n_pose @ Lxt 
    Rotation = _bend_constraint_b_get_Rt(S)
    #[a b c d e f g h i ......      a b c d e f g h i] shape : [Nx3x3, 1]
    new_Lxt = precomputed_rhs @ Rotation.reshape(-1,1)
    # new_Lxt = _bend_constraint_b_dense_part(Rotation, Lxt)
    b[...]=  new_Lxt.reshape(-1,1)
    


def stretching_constraint_A(ks : float ,v_size : int, edges ):
    row = []
    col = []
    data = []
        
    def append_block(row_offset, col_offset,x):
        nonlocal row, col, data
        row.append(row_offset+0);row.append(row_offset+1);row.append(row_offset+2)
            
        col.append(col_offset+0);col.append(col_offset+1);col.append(col_offset+2)
        data.append(x);data.append(x);data.append(x)

    append_data = make_sparse_matrix_triplet_function(row, col, data)
    row_size = len(edges)*3*2 # edge_size*pair*dimension
    # row_size = 3*v_size
    col_size = 3 * v_size
    
    for idx, (i, j) in enumerate(edges):
        # append_block( 3*i  , 3*i, 0.5 ); append_block(3*i , 3*j, -0.5 )
        # append_block( 3*j, 3*j, 0.5 ); append_block(3*j, 3*i, -0.5 )
        append_block( 6*idx + 3*0 , 3*i, 1.0 ); append_block(6 * idx + 3*0 , 3*j, -1.0 )
        append_block( 6*idx+3*1 , 3*j, 1.0 ); append_block(6 *idx + 3*1, 3*i, -1.0 )
    
    return sp.coo_matrix((data, (row, col)), shape = (row_size, col_size) , dtype=np.float64).tocsc()

def displacement_constraints_A(k_p : float, v_size : int):
    row = []
    col = []
    data = []
    
    append_data = make_sparse_matrix_triplet_function(row, col, data)
    col_size = v_size*3
    row_size = v_size*3
    for i in range(v_size):
        append_data(i, i, 1.0)
        # append_data(i, i, k_p)

    return sp.csc_matrix((data, (row, col)), shape = (row_size, col_size) , dtype=np.float64)


@njit
# def stretching_constraint_b(b,  v : np.ndarray, e, neutral_pose, neutral_rest_length : np.ndarray, __m_ks : float):
def stretching_constraint_b(b,  v : np.ndarray, e, neutral_pose, neutral_rest_length : np.ndarray):
        """
            e : edge index
            v : current mesh 


        """
        # e_idx1 = e[:, 0]
        # e_idx2 = e[:, 1]
        # v1 = v[e_idx1, : ]
        # v2 = v[e_idx2, : ]
        # spring = v2-v1 
        # edge_length = np.linalg.norm(spring,axis=-1 )
        # # spring_length = np.linalg.norm(spring)
        # # normalized_spring = spring / edge_length[..., None]
        # normalized_spring = spring / edge_length[..., None ]
        # delta =  (self.__m_nuetral_rest_stretch - edge_length ) * 0.5
        # direction = delta[..., None]*normalized_spring
        # pi = v1 + direction
        # pj = v2 - direction
        
        # for i, j in e:
            # b[3*i : 3*i+3, :] += self.__m_ks * 0.5 * ( pi[i, :] - pj[j, :] ).reshape(-1,1)
            # b[3*j : 3*j+3, :] += self.__m_ks * 0.5 * ( pj[j, :] - pi[i, :] ).reshape(-1,1)

        nv = neutral_pose
        for idx, (i, j) in enumerate(e):
            v1 = v[i, :]
            v2 = v[j, :]
            spring = v2 - v1
            # length = np.linalg.norm(spring, axis= - 1 )
            length = np.sqrt(np.sum(spring**2))

            normalized_spring = spring / length
            # n_length = np.linalg.norm(nv[j, :] - nv[i, :], axis=-1)
            n_length = np.sqrt(np.sum((nv[j, :] - nv[i, :])**2))
            delta = (length - n_length)*0.5
            
            
            
            
            direction = delta * normalized_spring
            pi = v1 + direction
            pj = v2 - direction
            # b[3*idx : 3*idx+3, :] +=   ( pi - pj ).reshape(-1,1)
            # b[6*idx+3*0 : 6*idx+3*0+3, :] +=   0.5*( pi - pj ).reshape(-1,1)
            # b[6*idx+3*1 : 6*idx+3*1+3, :] +=   0.5*( pj - pi ).reshape(-1,1)
            #orig
            b[6*idx+3*0 : 6*idx+3*0+3, :] +=   ( pi - pj ).reshape(-1,1)
            b[6*idx+3*1 : 6*idx+3*1+3, :] +=   ( pj - pi ).reshape(-1,1)
            
            # b[3*i : 3*i+3, :] +=  0.5* ( pi - pj ).reshape(-1,1)
            # b[3*j : 3*j+3, :] +=  0.5* ( pj - pi ).reshape(-1,1)




        
        

def displacement_constraints_b(b : sp.coo_matrix, k_p : sp.coo_matrix , neutral_v : np.ndarray, v : np.ndarray):
    """
        b = output
    """
    # assert(len(v) == len(neutral_v) and "len() size between neutral v and v is diff")


    # b += k_p@neutral_v.reshape(-1,1)
    b += neutral_v.reshape(-1,1)


def linearize_force( b,x_t, kp, ks, kb, neutral_mesh, edge , rest_stretch, bend_A, stretch_A, disp_A, precomuted_bend_const_rhs, precomputed_neibor, contact_forces_A, contact_tau_array, contact_coeff = 50.0):
    
    stretch_b = np.zeros((stretch_A.shape[0],1) )
    disp_b = np.zeros_like((len(neutral_mesh.v), 1))
    bend_b = np.zeros((bend_A.shape[0], 1))
    displacement_constraints_b(disp_b, kp, neutral_mesh.v, x_t)
    # stretching_constraint_b(stretch_b, x_t, edge, neutral_mesh.v , rest_stretch, ks)
    stretching_constraint_b(stretch_b, x_t, edge, neutral_mesh.v , rest_stretch)
    
    # bend_constraint_b(bend_b, x_t, bend_A, neutral_mesh, kb,precomuted_bend_const_rhs,precomputed_neibor)
    bend_constraint_b(bend_b, x_t, bend_A, neutral_mesh,precomuted_bend_const_rhs,precomputed_neibor)
    disp_b = disp_A.T @ kp @ disp_b 
    # stretch_b = stretch_A.T @ ks @ stretch_b
    stretch_b = stretch_A.T@ ks @ stretch_b
    # stretch_b = stretch_A.T@  stretch_b
    bend_b = bend_A.T@kb@bend_b
    contact_forces = contact_coeff*contact_forces_A.T@contact_tau_array
    # contact_forces = 1.0*contact_forces_A.T@contact_tau_array
    # b[...] = (disp_b + stretch_b + bend_b)
    b[...] = (disp_b + stretch_b + bend_b + contact_forces)
    # b[...] = (disp_b + stretch_b + bend_b )
    
def linearize_force_nonlinear( b, x_t, kp, ks, kb, neutral_mesh, edge , rest_stretch, bend_A, stretch_A, disp_A, precomuted_bend_const_rhs, precomputed_neibor, contact_forces_A, contact_tau_array, contact_coeff = 50.0):
    
    stretch_b = np.zeros((stretch_A.shape[0],1) )
    disp_b = np.zeros(( neutral_mesh.v.size,1))
    bend_b = np.zeros((bend_A.shape[0], 1))
    displacement_constraints_b(disp_b, kp, neutral_mesh.v, x_t)
    # stretching_constraint_b(stretch_b, x_t, edge, neutral_mesh.v , rest_stretch, ks)
    stretching_constraint_b(stretch_b, x_t, edge, neutral_mesh.v , rest_stretch)
    
    # bend_constraint_b(bend_b, x_t, bend_A, neutral_mesh, kb,precomuted_bend_const_rhs,precomputed_neibor)
    bend_constraint_b(bend_b, x_t, bend_A, neutral_mesh,precomuted_bend_const_rhs,precomputed_neibor)
    disp_b = disp_A.T @ kp @ disp_b 
    # stretch_b = stretch_A.T @ ks @ stretch_b
    stretch_b = stretch_A.T@ ks @ stretch_b
    # stretch_b = stretch_A.T@  stretch_b
    bend_b = bend_A.T@kb@bend_b
    contact_forces = contact_coeff*contact_forces_A.T@contact_tau_array
    # contact_forces = 1.0*contact_forces_A.T@contact_tau_array
    # b[...] = (disp_b + stretch_b + bend_b)
    b[...] = (disp_b + stretch_b + bend_b + contact_forces)
    # b[...] = (disp_b + stretch_b + bend_b )
    
    
@njit
def simulate_time_step(expression_pose, x_prev, u_t, phi, yt, step_size):
         
        exp = expression_pose
        x_t = phi @ u_t + yt
        x_acc_t_test = (x_t - x_prev.reshape(-1,1)) / step_size
        return x_t, x_acc_t_test
    
@njit
def simulate_time_step_nonlinear(x_prev, x_t, step_size):
        
        x_acc_t_test = (x_t.reshape(-1,1) - x_prev.reshape(-1,1)) / step_size
        return x_t, x_acc_t_test

# njit can't use sparse mat
def solve_phi_yt(bsums, B, sp_mass_matrix_inv, precomputed_Asums, prev_x_t_1, prev_x_acc, step_size):
    h = step_size
    h2 = step_size**2
    M_inv = sp_mass_matrix_inv
    damping_coeff_alpha = 0.99
    phi = precomputed_Asums(h2 * M_inv @ B)
    yt = precomputed_Asums(prev_x_t_1.reshape(-1,1) + damping_coeff_alpha*h*prev_x_acc.reshape(-1,1) + (h2*M_inv@bsums).reshape(-1,1) )
    return phi, yt
# njit can't use sparse mat

def solve_phi_yt_nonlinear(P, G, bsums, B, sp_mass_matrix_inv, precomputed_Asums, prev_x_t_1, prev_x_acc, step_size):
    h = step_size
    h2 = step_size**2
    M_inv = sp_mass_matrix_inv
    damping_coeff_alpha = 0.99
    phi = precomputed_Asums(h2 * P.T@G @ M_inv @ B )
    yy = (prev_x_t_1.reshape(-1,1)) + (damping_coeff_alpha*h*prev_x_acc.reshape(-1,1)) + (h2*M_inv@bsums).reshape(-1,1)
    yt = precomputed_Asums(P.T @G@ yy)
    return phi, yt
def solve_phi_yt_nonlinear2(P, G, bsums, B, sp_mass_matrix_inv, precomputed_Asums, prev_x_t_1, prev_x_acc, step_size):
    h = step_size
    h2 = step_size**2
    M_inv = sp_mass_matrix_inv
    damping_coeff_alpha = 0.99
    phi = precomputed_Asums(h2 * P.T@G @ M_inv @ B )
    yy = (prev_x_t_1.reshape(-1,1)) + (damping_coeff_alpha*h*prev_x_acc.reshape(-1,1)) + (h2*M_inv@bsums).reshape(-1,1)
    yt = precomputed_Asums(P.T @G@ yy)
    return phi, yt

def solve_ut( S, phi, yt, dt):
    

    S_Phi = S @ phi
    # S_Phi_T_S_Phi = S_Phi.T @ S_Phi

    dtSy = dt.reshape(-1,1) - S@yt
    # S_Phi_dtSy = S_Phi.T@(dtSy)

    result_u =np.linalg.lstsq(S_Phi, dtSy)[0]
    # result_u = np.linalg.solve(S_Phi_T_S_Phi, S_Phi_dtSy)
    return result_u

def solve_ut_nonlinear( reducedGxt, phi, yt, dt):
    

    # S_Phi = S @ phi
    S_Phi = phi
    
    # dtSy = dt.reshape(-1,1) - S@yt
    dtSy = reducedGxt + dt.reshape(-1,1) - yt

    result_u =np.linalg.lstsq(S_Phi, dtSy)[0]
    # result_u = np.linalg.solve(S_Phi_T_S_Phi, S_Phi_dtSy)
    return result_u


@njit
def static_solve(neutral, sel_exprs, marker_pose):
    A = sel_exprs
    neutral = neutral
    
    
    b = marker_pose.reshape(-1,1) - neutral
    

    w = np.linalg.solve( A.T@A, A.T @ b )
    w = np.clip(w, a_min=0.0, a_max=1.0)
    return w
            