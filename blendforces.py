import numpy as np 
import mesh as mm
import igl 
import typing 
import blendshapes
import scipy.sparse as  sp 
import geometry_helper as geo
from scipy.sparse.linalg import inv
from sksparse.cholmod import cholesky as spchol
import itertools
import collision_handler as ch
import logging, sys 
import gmath as gm 
import functools
from numba import njit
stiffnesses_args = typing.Union[tuple, str]
mass_args = typing.Union[typing.List[float], int, np.ndarray]
from scipy.sparse import spdiags

logger = logging.getLogger()

stream_handler = logging.StreamHandler(sys.stdout)
logger.addHandler(stream_handler)
logging.basicConfig(level=logging.DEBUG)

def make_sparse_matrix_triplet_function(raw_list, col_list, data_list):
    def append_data(i, j, val):
        raw_list.append(3*i + 0); col_list.append(3*j + 0) ; data_list.append(val)
        raw_list.append(3*i + 1); col_list.append(3*j + 1) ; data_list.append(val)
        raw_list.append(3*i + 2); col_list.append(3*j + 2) ; data_list.append(val)
    return append_data



class BlendForces:


    STRETCH_IDX :int
    BEND_IDX:int 
    DISP_IDX :int 

    def __init__(self, stiffnesses  : stiffnesses_args = "auto", iteration_num : int = 10, step_size = 0.16, mass : mass_args = 100.0, tau = 0.01):
        self.__m_stiffnesses = stiffnesses
        self.set_iteration_num(iteration_num)
        self.set_step_size(step_size)
        self.reset_simulation()
        self.__m_mass = mass

        self.__m_spatial_hash_object : ch.SpatialHashing = ch.SpatialHashing()

        self.__m_tau = 0.01

        data = np.load("D:\\lab\\2022\\mycode\\projective-dynamics-blendforce\\kpkskb.npz")
        dd = data['x'].reshape(-1)
        kp_data = dd[:6706]
        self.__m_kp = sp.kron(np.eye(3,3), spdiags(kp_data, 0, kp_data.size, kp_data.size))
        ks_data = dd[6706:26534]
        self.__m_ks = sp.kron(np.eye(3,3), spdiags(np.repeat(ks_data, 2), 0, ks_data.size*2, ks_data.size*2))
        kb_data = dd[26534:]
        self.__m_kb = sp.kron(np.eye(3,3), spdiags(kb_data, 0 , kb_data.size, kb_data.size))
        
        #TODO for testing
        # self.__m_kp = 1.0
        # self.__m_kp = 1.0
        # self.__m_ks = 1.0
        # self.__m_kb = 1.0


    def set_blendshapes(self, bld : blendshapes.Blendshapes):
        self.__m_blendshapes = bld

    @property
    def mesh(self):
        return self.__m_mesh
    
    @mesh.setter
    def mesh(self, mmm : mm.Mesh ):
        self.__m_mesh = mmm 
    

    def __compute_mass_diagonal(self):
        """
            make M matrix
        """
        pass 
    

    def __bend_force(self, v, sp_laplacian, neutral_mesh: mm.Mesh):
        Ai = self.__m_precomputed_constraint_Ai[BlendForces.BEND_IDX]
        
        b = np.zeros_like(self.__m_b)
        self.__bend_constraint_b(b, sp_laplacian, v, neutral_mesh)
        qq = Ai@v.reshape(-1,1)
        return (Ai@v.reshape(-1,1) - b)


    def __displace_force(self, v, neutral):
        return v.reshape(-1,1) - neutral.reshape(-1,1)
    
    def __stretch_force(self, v, neutral : mm.Mesh):
        e = neutral.e 
        v0 = v[e[:, 0], :]
        v1 = v[e[:, 1], :]
        
        
        vvv = (v0 - v1) - (neutral.v[e[:, 0], :] - neutral.v[e[:, 1], :])
        # return self.__stretch(v, neutral.e) - self.__stretch(neutral.v, neutral.e)
        return vvv.reshape(-1,1)


    def __stretch(self, v, e):
        """
            v : N x 3
            f : K x 3

            return 
                shapes : N
                [ e_length ] 
            
        """
        v0 = v[e[:, 0], :]
        v1 = v[e[:, 1], :]
        e1 = np.sqrt(np.sum((v0 - v1)**2, axis=-1))
        return e1



    
    def __displacement_constraints_b(self, b : np.ndarray, k_p : float , neutral_v : np.ndarray, v : np.ndarray):
        """
            b = output
        """
        assert(len(v) == len(neutral_v) and "len() size between neutral v and v is diff")


        b += k_p*neutral_v.reshape(-1,1)


    def __displacement_constraints_A(self, data : list,  row_indices : list, col_indices : list, k_p : float, neutral_v : np.ndarray):

        """
            # x_i : variable what we want to calculate
            # x_resp : rest pose(it means now, neutral pose of blendshapes)
            # W(x) = (x - x_resp)^2 
            # its derivative formula is W'(x) = 2*(x - x_ref) 
            # we want decarease its energy. minimum points naively 2*(x - x_ref) = 0 :=> x = x_ref
            # Its term be like here. 
            #   I x = b
        """
        append_data = make_sparse_matrix_triplet_function(row_indices, col_indices, data)

        for i in range(len(neutral.v)):
            append_data(i, i, k_p)

        return data, row_indices, col_indices

            
        

    def __bend_constraint_b(self, b : np.ndarray, sp_laplacian : sp.csc_matrix, vv : np.ndarray, neutral_mesh : mm.Mesh ):
        """

            vv : current vertices
            netural_vv : rest pose shape.
        """
        
        neutral_vv = neutral_mesh.v
        L = sp_laplacian

        

        S_i = np.empty((3,3), dtype=np.float64)
        self.RR = []
        for vidx_i in range(len(neutral_vv)):
            S_i[...] = 0.0
            for vidx_j in neutral_mesh.halfedge.v2v(vidx_i):
                wij = sp_laplacian[vidx_i,vidx_j]
                e_ij = vv[vidx_i] - vv[vidx_j]
                e_prime_ij = neutral_vv[vidx_i] - neutral_vv[vidx_j]
                S_i += wij*e_ij.reshape(-1,1) @ e_prime_ij.reshape(1,-1)
            u, _, vh = np.linalg.svd(S_i)
            R_i = vh.T@u
            # if np.linalg.det(R_i) == 0 :
                # u[-1, -1]
                # R_i = 

            self.RR.append(R_i)
        
        # reuse it 
        Rij = S_i

        index = 0 
        b_coeff = np.empty((3,1))
        for vidx_i in range(len(neutral_mesh.v)):
            Rij[...]  = 0.0
            b_coeff[...] = 0.0
            for vidx_j in neutral_mesh.halfedge.v2v(vidx_i):
                w_ij = -sp_laplacian[vidx_i, vidx_j ] #
                Rij = (self.RR[vidx_i] + self.RR[vidx_j])
                e_ij = vv[vidx_i] - vv[vidx_j]
            
                b_coeff += w_ij*(Rij @ e_ij.reshape(-1,1))
            b[index*3:index*3+3, :]  = b_coeff[...]

    def __bend_constraint_A(self, data : list, row_indices : list, col_indices : list, sp_laplacian : sp.csc_matrix, neutral_vv ):
        """
        
            orga sorkine 
            https://igl.ethz.ch/projects/ARAP/arap_web.pdf
            bending constaints
        """
        row, col = sp_laplacian.nonzero()
        dt = sp_laplacian.data

        append_data = make_sparse_matrix_triplet_function(row_indices, col_indices, data)
        # row_offset =row_indices[-1] + 1
        # row += row_offset 
        data.extend(dt)
        row_indices.extend(row)
        col_indices.extend(col)
        return row_indices, col_indices, data
        
    def __stretching_constraint_b(self, b, v : np.ndarray, e, neutral_rest_length : np.ndarray, __m_ks : float):
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

        nv = self.__m_blendshapes.neutral_pose()
        for i, j in e:
            v1 = v[i, :]
            v2 = v[j, :]
            spring = v2 - v1
            length = np.linalg.norm(spring, axis= - 1 )
            normalized_spring = spring / length
            n_length = np.linalg.norm(nv[j, :] - nv[i, :], axis=-1)
            delta = (n_length - length)*0.5
            direction = delta[..., None] * normalized_spring
            pi = v1 + direction
            pj = v2 - direction
            b[3*i : 3*i+3, :] += __m_ks * 0.5 * ( pi - pj ).reshape(-1,1)
            b[3*j : 3*j+3, :] += __m_ks * 0.5 * ( pj - pi ).reshape(-1,1)




        
        

    def __stretching_constraint_A(self, data :list, row_indices : list, col_indices : list, ks : float ,v : np.ndarray, edges ):
        """
            ocasionally author called it. elastic forces or Stretch.
            W_stretch(x) we simply shorten it as Ws(x)
            Ws(x) = Sum { (||x_i - x_j || - r_ij)^2 } 
            derivate it by edges e_ij := x_i - x_j
            then we can get it derivative Ws'(x)
                Ws'(x) = Sum { x_i - x_j } - Sum { (r_ij*(x_i-x_j))/||x_i - x_j|| }

            find its solution on zero gradient position. Sum { x_i - x_j } - Sum { (r_ij*(x_i-x_j))/||x_i - x_j= || }
            
            Do not use this term solely, It may occur sigular matrix problem.(It occured when I tested it by using three points of triangle.)

            It looks simillar with L(Laplacian matrix.) 
            But Laplacian matrix was not made for and used for considering length between vertice(called Edge).

            Laplacian used for blending and smooth curvature.

        """
        
        append_data = make_sparse_matrix_triplet_function(row_indices, col_indices, data)



        # for i, j in edges:
            
        #     # append_data(i, i, -1.0 * self.__m_ks * 0.5); append_data(i, j, 1.0 * self.__m_ks * 0.5)
        #     # append_data(j, j, -1.0 * self.__m_ks * 0.5); append_data(j, i, 1.0 * self.__m_ks * 0.5)
            
        #     # before below 
        #     # append_data(i, i, -1.0 * self.__m_ks* 0.5); append_data(i, j, 1.0 * self.__m_ks* 0.5)
        #     # append_data(j, j, -1.0 * self.__m_ks* 0.5); append_data(j, i, 1.0 * self.__m_ks* 0.5)
            
        #     append_data(i, i, 1.0 * self.__m_ks* 0.5); append_data(i, j, -1.0 * self.__m_ks* 0.5)
        #     append_data(j, j, 1.0 * self.__m_ks* 0.5); append_data(j, i, -1.0 * self.__m_ks* 0.5)
            
            
        for idx ,(i, j) in enumerate(edges):
            append_data(2*idx + 0, i ,1.0 * self.__m_ks), append_data(2*idx + 0, j,-1.0 * self.__m_ks)
            append_data(2*idx + 1, j, 1.0 * self.__m_ks), append_data(2*idx + 1, i, -1.0 * self.__m_ks)

    

    def __contact_reponse_b(self):
        pass
    def __contact_reponse_A(self, data, row, col):
        v_indice = self.__m_spatial_hash_object.query()

        



    

    def __compute_mass_matrix(self):
        N, _ = self.__m_blendshapes.neutral_pose().shape

        diag = np.empty((N*3), dtype=np.float64)
        diag[...] = self.__m_mass
        self.__m_sp_mass_matrix =  sp.spdiags(diag, 0, diag.size, diag.size).tocsc()
        self.__m_sp_mass_matrix_inv = sp.spdiags(1.0/diag, 0, diag.size, diag.size).tocsc()
        return self.__m_sp_mass_matrix

    def __precompute_neutral_pose_constants(self):
        self.__compute_mass_matrix()
        self.__m_nuetral_rest_stretch = self.__stretch(self.__m_blendshapes.neutral_pose(), self.__m_blendshapes.neutral_mesh().e)
        self.__m_sp_laplacian_matrix = geo.make_laplacian(self.__m_blendshapes.neutral_mesh())


        self.__m_bs_expression_matrix = self.__m_blendshapes.expression_pose()
        self.__m_bs_netural_pose = self.__m_blendshapes.neutral_pose().reshape(-1,1)

        self.__m_bs_selectec_marker_expression_matrix = self.__m_blendshapes.expression_pose(self.__m_marker_indices)
        self.__m_bs_selected_netural_pose = self.__m_blendshapes.neutral_pose()[self.__m_marker_indices, :].reshape(-1,1)



    def set_step_size(self, size : float):
        self.__m_step_size = size

    def set_damping_factor(self, damp_factor :float):
        self.__m_damping_factor = max(min(damp_factor, 1.0), 0.0)

        

        


    def set_iteration_num(self, num : int):
        if num < 0 :
            raise ValueError("iteration must be signed integer")
        self.__m_iteration_num  = num 

    def reset_simulation(self):
        self.__m_first_iter_flag  = True 
        self.__m_As_sum = None 
        self.__m_bs_sum = None



    def precompute(self):
        self.__precompute_neutral_pose_constants()

        self.__m_vN, *_ = self.__m_blendshapes.neutral_pose().shape
        self.__m_b = np.empty((3*self.__m_vN,1), dtype=np.float64)
        

        self.__m_Ai = sp.identity(self.__m_vN, dtype=np.float64)
        self.__m_Bi = self.__m_Ai  # Ai == Bi (only affect numerical solution procedure.)
        print(self.__m_blendshapes.neutral_mesh().e)
        func_list = []
        #TODO
        
        BlendForces.STRETCH_IDX = 0
        BlendForces.BEND_IDX = 1
        BlendForces.DISP_IDX = 2
        self.__m_w_coeff, self.__m_sp_neighbor_mat = geo.make_neibor_matrix(self.__m_blendshapes.neutral_mesh())
        vvv = self.__m_sp_neighbor_mat@self.__m_blendshapes.neutral_mesh().v
        self.__m_weighted_neutral_edges = self.__m_w_coeff.reshape(-1,1)* vvv
        self.__m_precomuted_bend_const_rhs = gm.precompute_arap_synbolic_S_sums(self.__m_blendshapes.neutral_mesh())
        self._m_sp_linearize_arap_mat = gm.linearize_arap_params(self.__m_blendshapes.neutral_mesh())
        
        self.__m_precomputed_constraint_Ai = {}
        self.__m_precomputed_constraint_Ai[BlendForces.BEND_IDX] = {"coeff" : self.__m_kb, "A" : gm.bend_constraint_A( self.__m_blendshapes.neutral_mesh(), self.__m_kb)}
        self.__m_precomputed_constraint_Ai[BlendForces.DISP_IDX] = {"coeff" : self.__m_kp, "A" : gm.displacement_constraints_A(-self.__m_kp, self.__m_vN)}
        self.__m_precomputed_constraint_Ai[BlendForces.STRETCH_IDX] = {"coeff" : self.__m_ks, "A" : gm.stretching_constraint_A(-self.__m_ks, self.__m_vN, self.__m_blendshapes.neutral_mesh().e)}

        self.__m_As_sum = sp.csc_matrix((self.__m_vN*3, self.__m_vN*3) , dtype=np.float64)
        for key, val in  self.__m_precomputed_constraint_Ai.items():
            if(val is not None ):
                coeff = val["coeff"]
                A = val["A"]
                self.__m_As_sum += -A.T@ coeff @ A
        
        # d_A = self.__m_precomputed_constraint_Ai[BlendForces.DISP_IDX]['A']
        # c_d =   self.__m_precomputed_constraint_Ai[BlendForces.DISP_IDX]['coeff']
        # s_A = self.__m_precomputed_constraint_Ai[BlendForces.STRETCH_IDX]['A']
        # c_s = self.__m_precomputed_constraint_Ai[BlendForces.STRETCH_IDX]['coeff']
        # b_A = self.__m_precomputed_constraint_Ai[BlendForces.BEND_IDX]['A']
        # c_b = self.__m_precomputed_constraint_Ai[BlendForces.BEND_IDX]['coeff']
        # self.__m_As_sum = -(d_A.T@c_d@d_A)
        # self.__m_As_sum = -(s_A.T@c_s@s_A)
        # self.__m_As_sum = -(b_A.T@c_b@b_A)

        I = sp.identity(self.__m_sp_mass_matrix.shape[0]).tocsc()
        h = self.__m_step_size
        h2 = self.__m_step_size**2
        M_inv = self.__m_sp_mass_matrix_inv
        
        tmp1 = self._tmp1 = (I -  h2 * M_inv@ self.__m_As_sum)
        
        # phi = tmp1 @ tmp2 
        self.__m_precomputed_I_Asums2 = tmp1
        self.__m_precomputed_I_Asums = spchol(tmp1)
        print("precompute")


    def add_marker_index(self, index_list : typing.Union[typing.List[int], np.ndarray]):
        self.__m_marker_indices = index_list
        # data = [1 for _ in (self.__m_marker_indices)]
        # row = [ i for i in range(len(self.__m_marker_indices))]
        # col = self.__m_marker_indices

        data = []
        rows = [] 
        cols = []
        
        marker_len = len(self.__m_marker_indices)
        # append_data = make_sparse_matrix_triplet_function(rows, cols, data)

        row = 0 
        def append_data(i, j, val):
            rows.append(3*i     ); cols.append(3*j + 0) ; data.append(val)
            rows.append(3*i + 1 ); cols.append(3*j + 1) ; data.append(val)
            rows.append(3*i + 2 ); cols.append(3*j + 2) ; data.append(val)
        for i in  self.__m_marker_indices:
            append_data(row    , i  , 1)
            row += 1
        
        N = len(self.__m_blendshapes.neutral_pose())
        self.__M_S_sp_mat = sp.csc_matrix( (data, (rows, cols)), shape=(marker_len*3, N*3), dtype=np.float64) 
        



    def __linearlize_forces(self, x_t):
        """
            see appendix
            A_i = -k*F_i.T@F_i
            b_i = -k*F_i.T@G_i@p_i  
        """
        #disp cons 
        self.__m_b[...] = 0
        
        # self.__m_bd = np.zeros_like(self.__m_b)
        self.__displacement_constraints_b(self.__m_b, self.__m_kp, self.__m_blendshapes.neutral_pose(), x_t)
        # self.__displacement_constraints_b(self.__m_bd, self.__m_kp, self.__m_blendshapes.neutral_pose(), x_t)
        # A1, b1 = self.__solve_projective_dynamaics(self.__m_precomputed_constraint_Ai[BlendForces.DISP_IDX], self.__m_b, self.__m_ks)
        # self.__m_b += self.m_As[1].transpose()@self.__m_bd
        #stretching cons 
        # self.__m_kb
        # self.__m_sb = np.zeros_like(self.__m_b)
        self.__stretching_constraint_b(self.__m_b, x_t, self.__m_blendshapes.neutral_mesh().e, self.__m_nuetral_rest_stretch, self.__m_ks)
        # self.__stretching_constraint_b(self.__m_sb , x_t, self.__m_blendshapes.neutral_mesh().e, self.__m_nuetral_rest_stretch, self.__m_ks)
        # self.__m_b += self.m_As[0].transpose()@self.__m_sb

        #bending cons 
        # self.__bend_constraint_b(self.__m_b, self.__m_sp_laplacian_matrix, x_t ,self.__m_blendshapes.neutral_mesh())
        


        # conatact(collision) cons  
        # TODO 


        
        
        return self.__m_b




    def __static_solve(self, marker_pose):
        A = self.__m_bs_selectec_marker_expression_matrix
        neutral = self.__m_bs_selected_netural_pose
        
        
        b = marker_pose.reshape(-1,1) - neutral
        

        w = np.linalg.solve( A.T@A, A.T @ b )
        w = np.clip(w, a_min=0.0, a_max=1.0)
        return self.__m_blendshapes.make_pose_by_weight(w)
        

    def solve_phi_yt(self, Asums, bsums, B, prev_x_t_1, prev_x_acc):
        """
            Asum : precomuted strecth, disp, bend, contact constraints coeff
            bsum : linearize forces
            B : blendshapes
            prev_xt_t_1 : prev pos
            prev_x_acc : prev acc
        """

        h = self.__m_step_size
        h2 = self.__m_step_size**2
        M_inv = self.__m_sp_mass_matrix_inv
        a = 0.1
        phi = self.__m_precomputed_I_Asums(h2 * M_inv @ B)
        yt = self.__m_precomputed_I_Asums(prev_x_t_1.reshape(-1,1) + a*h*prev_x_acc.reshape(-1,1) + (h2*M_inv@bsums).reshape(-1,1) )
        return phi, yt


    def solve_ut(self, phi, yt, dt):
        S = self.__M_S_sp_mat

        S_Phi = S @ phi
        S_Phi_T_S_Phi = S_Phi.T @ S_Phi

        dtSy = dt.reshape(-1,1) - S@yt
        S_Phi_dtSy = S_Phi.T@(dtSy)

        result_u =np.linalg.lstsq(S_Phi, dtSy)[0]
        # result_u = np.linalg.solve(S_Phi_T_S_Phi, S_Phi_dtSy)
        return result_u
        

    
        
        


    


    def simulate_time_step(self, x_prev, x_acc_prev, u_t, phi, yt):
         
        exp = self.__m_blendshapes.expression_pose()
        # x_acc_t = x_acc_prev + self.__m_step_size * (self.__m_sp_mass_matrix_inv @ (exp @ u_t)).reshape(-1, 3)
        # x_t = x_prev + self.__m_step_size*x_acc_t
        ###

        x_t = phi @ u_t + yt
        # x_acc_t = x_acc_prev.reshape(-1,1) + self.__m_step_size * self.__m_sp_mass_matrix_inv @ (exp@u_t + self.__m_As_sum @ x_t + self.__m_b)
        x_acc_t_test = (x_t - x_prev.reshape(-1,1)) / self.__m_step_size
        return x_t, x_acc_t_test
        
        
        # self.__m_As_sum @ appx_x_t + 

        # force = (fext + sys_f)
        # s = (self.__m_step_size*(self.__m_sp_mass_matrix_inv ) @ force)
        # print(s.shape)
        # x_acc_t = x_acc_prev + (s).reshape(-1,3)
        # x_t = x_prev + self.__m_step_size * x_acc_t
        return x_t, x_acc_t
    
    def update2(self, new_marker_pos, frame, contact_A : sp.csc_matrix, contact_tau_array :np.ndarray, contact_coeff = 50.0):
        
        if self.__m_first_iter_flag:
            self.__m_first_iter_flag = False 
            w  = gm.static_solve(self.__m_bs_selected_netural_pose, self.__m_bs_selectec_marker_expression_matrix,  new_marker_pos)
            self.__x_prev = self.__m_blendshapes.make_pose_by_weight(w)
            self.__x_acc_prev = np.zeros_like(self.__x_prev)
        
        M_inv = self.__m_sp_mass_matrix_inv
        I = sp.identity(self.__m_sp_mass_matrix.shape[0]).tocsc()
        h2 = self.__m_step_size**2
        M_inv = self.__m_sp_mass_matrix_inv
        tmp1 = self._tmp1 = (I -  h2 * M_inv@ (self.__m_As_sum  +(-contact_coeff*contact_A.T@contact_A)))
        # tmp1 = self._tmp1 = (I -  h2 * M_inv @ (self.__m_As_sum))
        self.__m_precomputed_I_Asums = spchol(tmp1)
        frames = []
        x_t = self.__x_prev + self.__m_step_size* self.__x_acc_prev
        for f in range(frame):
            for iter_n in range(self.__m_iteration_num):
                # logger.debug("%d", iter_n)
                gm.linearize_force(self.__m_b, x_t, self.__m_kp, self.__m_ks,self.__m_kb, \
                    self.__m_blendshapes.neutral_mesh(), \
                        self.__m_blendshapes.neutral_mesh().e,self.__m_nuetral_rest_stretch, \
                            self.__m_precomputed_constraint_Ai[BlendForces.BEND_IDX]["A"], \
                                self.__m_precomputed_constraint_Ai[BlendForces.STRETCH_IDX]["A"], \
                                    self.__m_precomputed_constraint_Ai[BlendForces.DISP_IDX]["A"] ,\
                                        self.__m_precomuted_bend_const_rhs,self._m_sp_linearize_arap_mat, contact_A, contact_tau_array) 
                
                self.phi, self.y_t =gm.solve_phi_yt(self.__m_b, self.__m_bs_expression_matrix, self.__m_sp_mass_matrix_inv, self.__m_precomputed_I_Asums, self.__x_prev, self.__x_acc_prev, self.__m_step_size)
                u_t = gm.solve_ut(self.__M_S_sp_mat, self.phi, self.y_t, new_marker_pos)
                x_t, x_acc_t = gm.simulate_time_step(self.__m_blendshapes.expression_pose(), \
                    self.__x_prev, u_t,  self.phi, self.y_t, self.__m_step_size)
                
                # self.phi, self.y_t = self.solve_phi_yt(self.__m_As_sum, self.__m_b, \
                                            # B = self.__m_bs_expression_matrix, \
                                            # prev_x_t_1= self.__x_prev, prev_x_acc= self.__x_acc_prev)
                # u_t = self.solve_ut(self.phi, self.y_t, new_marker_pos)
                # x_t, x_acc_t = self.simulate_time_step(self.__x_prev, self.__x_acc_prev, u_t,  self.phi, self.y_t)
                x_t, x_acc_t = x_t.reshape(-1,3), x_acc_t.reshape(-1,3)
                self.__x_prev, self.__x_acc_prev = x_t, x_acc_t
            frames.append(x_t)
        return frames
    
    def update(self,  new_marker_pos, frame):
        if self.__m_first_iter_flag:
            self.__m_first_iter_flag = False 
            self.__x_prev = self.__static_solve(new_marker_pos)
            self.__x_acc_prev = np.zeros_like(self.__x_prev)
        # self.__m_first_iter_flag = False 
        # self.__x_prev = self.__static_solve(new_marker_pos)
        # self.__x_acc_prev = np.zeros_like(self.__x_prev)

        frames = []
        x_t = self.__x_prev + self.__m_step_size* self.__x_acc_prev
        for f in range(frame):
            for iter_n in range(self.__m_iteration_num):
                logger.debug("%d", iter_n)
                self.__linearlize_forces(x_t) # update self.__m_b 
                
                self.phi, self.y_t = self.solve_phi_yt(self.__m_As_sum, self.__m_b, \
                                            B = self.__m_bs_expression_matrix, \
                                            prev_x_t_1= self.__x_prev, prev_x_acc= self.__x_acc_prev)
                u_t = self.solve_ut(self.phi, self.y_t, new_marker_pos)
                # u_t = np.clip(u_t, 0.0, 1.0)
                # print(u_t)
                x_t, x_acc_t = self.simulate_time_step(self.__x_prev, self.__x_acc_prev, u_t,  self.phi, self.y_t)
                x_t, x_acc_t = x_t.reshape(-1,3), x_acc_t.reshape(-1,3)
                
                # x_t = self.__m_blendshapes.make_pose_by_weight(u_t).reshape(-1,3)
                # x_acc_t = np.zeros_like(self.__x_prev) 
                
            self.__x_prev, self.__x_acc_prev = x_t, x_acc_t
            frames.append(x_t)
        # return x_t
        return frames


if __name__ == "__main__":
    import os , glob 
    import cv2
    import subprocess
    import numpy as np
    
    # MediaPipe Face Mesh 초기화
    # mp_face_mesh = mp.solutions.face_mesh
    # face_mesh = mp_face_mesh.FaceMesh(static_image_mode=False, max_num_faces=1, refine_landmarks=True)

    # dlib 68개 landmark 인덱스 (MediaPipe 468개 포인트 중 dlib에 해당하는 것만 선택)
    DLIB_68_IDX = [162,234,93,58,172,136,149,148,152,377,378,365,397,288,323,454,389,71,63,105,66,107,336,
                    296,334,293,301,168,197,5,4,75,97,2,326,305,33,160,158,133,153,144,362,385,387,263,373,
                    380,61,39,37,0,267,269,291,405,314,17,84,181,78,82,13,312,308,317,14,87]


    
    data_path = "D:\\lab\\2022\\mycode\\FaceCaptureWithIK\\data\\ICT-data"
    neutral_pth = os.path.join(data_path, "generic_neutral_mesh.obj")
    shapes_path = os.path.join(data_path, "shapes")
    file_pths = glob.glob(os.path.join(shapes_path, "**.obj"))
    neutral = mm.Mesh()
    neutral.load_from_file(neutral_pth)
    bs_list = []
    for fpth in file_pths:
        m = mm.Mesh()
        m.load_from_file(fpth)
        bs_list.append(m)
        
    import asyncio
    
    lmk_idx = [1278,1272,12,1834,243,781,2199,1447,966,3661,4390,3022,2484,4036,2253,3490,3496,268,493,1914,2044,1401,3615,4240,4114,2734,2509,978,4527,4942,4857,1140,2075,1147,4269,3360,1507,1542,1537,1528,1518,1511,3742,3751,3756,3721,3725,3732,5708,5695,2081,0,4275,6200,6213,6346,6461,5518,5957,5841,5702,5711,5533,6216,6207,6470,5517,5966,]
    
    datas = []
    for pth in glob.glob("./exported_objs/**.obj"):
        m = mm.Mesh()
        m.load_from_file(pth)
        datas.append(m.v[lmk_idx, :])
        # datas.append(m.v)
    
    # sub = subprocess.Popen("python ./viewer.py", stdin = subprocess.PIPE)


    bshapes = blendshapes.Blendshapes(neutral, bs_list)
    bshapes.build()
    bb = BlendForces()
    bb.set_blendshapes(bld=bshapes)
    bb.add_marker_index(lmk_idx)
    bb.precompute()
    xt = neutral[lmk_idx]
    # for it in range(50) : 
    a = (np.sin(np.linspace(0, 100, 3000))).reshape(-1,1)
    av = np.hstack([a, np.zeros_like(a), np.zeros_like(a)])
    ii = 0 


    def __static_solve(marker_pose):
        A = bshapes.expression_pose(lmk_idx)
        neutral = bshapes.neutral_pose()[lmk_idx, :]
        
        
        b = marker_pose.reshape(-1,1) - neutral.reshape(-1,1)
        
        w = np.linalg.solve( A.T@A, A.T @ b )
        # w = np.linalg.lstsq(A,b)[0]
        w = np.clip(w, a_min=0.0, a_max=1.0)
        bb = bshapes.make_pose_by_weight(w)
        return bb
    w = np.zeros((len(bs_list), 1))
    print(len(bs_list))
    
    
    
    async def run_sub(queue, loop)    :
        print(loop, "ewew")
        proc = await loop.run_in_executor(
                None, 
                subprocess.Popen, 
                "python ./viewer.py", 0, None,
                subprocess.PIPE,       # stdin,
                
            )
        # proc = subprocess.Popen("python ./viewer.py", stdin=subprocess.PIPE)
        ii = 0 
        
        while True:
            data = await queue.get()
            proc.stdin.write(data)
            proc.stdin.flush()
            # print(f"put {ii}")
            ii+=1


    
    def run_main(queue, loop):
        # cap = cv2.VideoCapture(0)
        import spatialhashing as sph
        marker = neutral[lmk_idx] 
        ii = 0
        aai = ii % len(datas)
        marker = datas[ii]
        import time
        hashgrid = sph.OptimSpatialHashGrid()
        hashgrid.attach_face_data(neutral.f)
        hashgrid.upate_vertex_data(neutral.v)
        frames = []
        datasize = len(datas)
        for i, marker in enumerate(datas)  : 
            print(f"{i} th/%{datasize}")
            candidates = hashgrid.query_overlapped_tris()
            e, contact_A , tau_array = hashgrid.calc_forces(candidates)
            frame = bb.update2(marker, frame=1, contact_A=contact_A, contact_tau_array=tau_array)[0]
            
            frame[lmk_idx] = marker
            frames.append(frame)
            hashgrid.upate_vertex_data(frame)
            asyncio.run_coroutine_threadsafe(queue.put(frame), loop)

        if not os.path.exists("./testanim_out"):
            os.makedirs("./testanim_out")
        np.save("./testanim_out/testanim.npy",  np.array(frames))
           
        marker = neutral[lmk_idx] 
        ii = 0
        aai = ii % len(datas)
        marker = datas[ii]
        import time
        while True : 
            start = time.time()
            candidates = hashgrid.query_overlapped_tris()
            f = hashgrid.calc_forces(candidates)
            frame = bb.update2(marker, frame=1)[0]

            end = time.time()
            
            print(f'{end - start} f/s')
            frame[lmk_idx] = datas[aai]
            ii += 1
            aai = ii % len(datas)
            marker = datas[aai]
       
            asyncio.run_coroutine_threadsafe(queue.put(frame), loop)
           

    import threading
    queue = asyncio.Queue()
    loop = asyncio.get_event_loop()
    print(loop, "init")
    calculation_thread = threading.Thread(target=run_main, args=(queue,loop))
    calculation_thread.daemon = True  # Ensure the thread exits when the main program exits
    calculation_thread.start()

    loop.run_until_complete(run_sub(queue,loop))
    
    
    import multiprocessing as mp 
    
    