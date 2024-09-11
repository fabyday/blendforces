import numpy as np 
from . import mesh as mm
import igl 
import typing 
import blendshapes
import scipy.sparse as  sp 
import geometry_helper as geo
from scipy.sparse.linalg import inv
from sksparse.cholmod import cholesky as spchol


import functools

stiffnesses_args = tuple | typing.Literal["auto"] 
mass_args = typing.Union[typing.List[float], int, np.ndarray]
class BlendForces:


    STRETCH_IDX = 0
    BEND_IDX = 1
    DISPLACE_IDX = 2

    def __init__(self, stiffnesses  : stiffnesses_args = "auto", iteration_num : int = 10, step_size = 0.16, mass : mass_args = 10.0, tau = 0.01):
        self.__m_stiffnesses = stiffnesses
        self.set_iteration_num(iteration_num)
        self.set_step_size(step_size)
        self.reset_simulation()
        self.__m_mass = mass


        self.__m_tau = 0.01


        #TODO for testing
        self.__m_kp = 1.0
        self.__m_ks = 1.0
        self.__m_kb = 1.0


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
    

    

    def __stretch(self, v, f):
        """
            v : N x 3
            f : K x 3

            e1 e2 e3 is column vectors
            return 
                shapes : N x 3
                [ e0_length | e1_length | e3_length ] 
            
        """
        v0 = v[f[:, 0], :]
        v1 = v[f[:, 1], :]
        v2 = v[f[:, 2], :]
        e1 = np.sqrt(np.sum((v0 - v1)**2, axis=-1))
        e2 = np.sqrt(np.sum( (v1 - v2) **2, axis=-1))
        e3 = np.sqrt(np.sum( (v2 - v0) **2, axis=-1))
        return np.hstack([e1, e2, e3])




    def __displacement_constraints_b(self, b : np.ndarray, k_p : float , neutral_v : np.ndarray, v : np.ndarray):
        """
            b = output
        """
        assert(len(v) == len(neutral_v), "len() size between neutral v and v is diff")


        b += neutral_v.shape(-1,1)


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
        new_ind_row = row_indices[-1] + 1
        data += [1.0 for _ in range(len(neutral_v))]
        row_indices += [new_ind_row + i for i in range(len(neutral_v))]
        col_indices += [ i for i in range(len(neutral_v))]
        

        return data, row_indices, col_indices

            
        

    def __bend_constraint_b(self, b : np.ndarray, sp_laplacian : sp.csc_matrix, vv : np.ndarray, neutral_mesh : mm.Mesh ):
        """

            vv : current vertices
            netural_vv : rest pose shape.
        """

        neutral_vv = neutral_vv.v
        S_i = np.empty((3,3), dtype=np.float64)
        RR = []
        for vidx_i in range(len(neutral_vv.v)):
            S_i[...] = 0.0
            for vidx_j in neutral_vv.halfedge.v2v(vidx_i):
                wij = sp_laplacian[vidx_i,vidx_j]
                e_ij = vv[vidx_i] - vv[vidx_j]
                e_prime_ij = neutral_vv[vidx_i] - neutral_vv[vidx_j]
                S_i += wij*e_ij.reshape(-1,1) @ e_prime_ij.reshape(1,-1)
            u, _, vh = np.linalg.svd(S_i)
            R_i = vh.T@u
            # if np.linalg.det(R_i) == 0 :
                # u[-1, -1]
                # R_i = 

            RR.append(R_i)
        
        # reuse it 
        Rij = S_i

        index = 0 
        b_coeff = np.empty((3,1))
        for vidx_i in range(len(neutral_vv.v)):
            Rij[...]  = 0.0
            b_coeff[...] = 0.0
            for vidx_j in neutral_vv.halfedge.v2v(vidx_i):
                w_ij = -sp_laplacian[vidx_i, vidx_j ] #
                Rij = (RR[vidx_i] + RR[vidx_j])
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


        row_offset =row_indices[-1] + 1
        row += row_offset 
        data.extend(dt)
        row_indices.extend(row)
        col_indices.extend(col)
        return 
        
    def __stretching_constraint_b(self, b, v : np.ndarray, f, neutral_rest_length : np.ndarray):
        edges = self.__stretch(v, f)
        
        b = neutral_rest_length/ edges
        

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
        # elastic forces.

        new_row_offset = row_indices[-1] + 1
        data += len(edges)*[ks*1,ks*(-1)]

        for i, j in edges:
            row_indices.append(new_row_offset)
            col_indices.append(i)
            col_indices.append(j)
            new_row_offset += 1


    

    def __contact_reponse_b(self):
        pass
    def __contact_reponse_A(self):
        pass



    

    def __compute_mass_matrix(self):
        N, _ = self.__m_blendshapes.neutral_pose().shape

        diag = np.empty((N), dtype=np.float64)
        diag[...] = self.__m_mass
        self.__m_sp_mass_matrix =  sp.sparse.spdiags(diag, 0, diag.size, diag.size)
        return self.__m_sp_mass_matrix

    def __precompute_neutral_pose_constants(self):
        self.__compute_mass_matrix()
        self.__m_nuetral_rest_stretch = self.__stretch(self.__m_mesh.v, self.__m_mesh.f)
        self.__m_sp_laplacian_matrix = geo.make_laplacian(self.__m_blendshapes.neutral_pose())


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
        self.__neutral_pose_rest_length, self.__lapalcian_matrix = self.__precompute_neutral_pose_constants()

        self.__m_vN, *_ = self.__m_blendshapes.neutral_pose().shape
        self.__m_b = np.empty((3*self.__m_vN,1), dtype=np.float64)
        

        self.__m_Ai = sp.identity(self.__m_vn, dtype=np.float64)
        self.__m_Bi = self.__m_Ai  # Ai == Bi (only affect numerical solution procedure.)
    
        func_list = []
        func_list += [functools.partial(self.__stretching_constraint_A, ks=self.__m_ks, v = self.__m_blendshapes.neutral_pose(), edges = )]
        func_list += [functools.partial(self.__bend_constraint_A,sp_laplacian = self.__m_sp_laplacian_matrix, neutral_vv  = self.__m_blendshapes.neutral_pose() )]
        func_list += [functools.partial(self.__displacement_constraints_A, k_p = self.__m_kp, neutral_v = self.__m_blendshapes.neutral_pose())]

        self.__m_precomputed_constraint_Ai = []
        for func in func_list:
            data = []
            row = []
            col = []
            func(data, row, col)
            Ai = sp.csc_matrix((data, row, col), shape =(self.__m_vN*3, self.__m_vn*3) , dtype=np.float64)
            self.__m_precomputed_constraint_Ai.append(Ai)
            



    def add_marker_index(self, index_list : list[int]):
        self.__m_marker_indices = index_list
        data = [1 for _ in len(self.__m_marker_indices)]
        row = [ i for i in range(len(self.__m_marker_indices))]
        col = self.__m_marker_indices

        N = len(self.__m_blendshapes.neutral_pose())
        self.__M_S_sp_mat = sp.csc_matrix((data, (row,col)), shape=(len(self.__m_marker_indices),N), dtype=np.float64) 
        



    def __solve_projective_dynamaics( A : sp.csc_matrix, b : np.ndarray, x : np.ndarray, k : float):
        """
            k : factor
        """
        
        chol = chol(-k*(A.T))
        pi = chol(-k*(A.T@A@x + b))

        Ai = -k*A.T@A 
        bi = -k*A.T@pi        

        return Ai, bi

    def __linearlize_forces(self, x_t):
        """
            see appendix
            A_i = -k*F_i.T@F_i
            b_i = -k*F_i.T@G_i@p_i  
        """
        data = [] 
        row_indices = []
        col_indices = []


        self.__m_b[...] = 0 # reset zero
        self.__displacement_constraints_b(data,  row_indices, col_indices, self.__m_b, self.__m_blendshapes.neutral_pose(), x_t)
        A1, b1 = self.__solve_projective_dynamaics(self.__m_Ai[0], self.__m_b)

        self.__m_b[...] = 0 # reset zero
        self.__stretching_constraint_b(data,  row_indices, col_indices, self.__m_b)
        A2, b2 = self.__solve_projective_dynamaics(self.__m_Ai[1], self.__m_b)

        self.__solve_sparse(self.__m_stiffnesses,)
         
        self.__m_b[...] = 0 # reset zero
        self.__bend_constraint_b(data,  row_indices, col_indices, self.__m_b)
        A3, b3 = self.__solve_projective_dynamaics(self.__m_Ai[1], self.__m_b)

        

        
        return A1 + A2 + A3, b1 + b2 + b3




    def __static_solve(self, marker_pose):
        A = self.__m_blendshapes.expression_pose()
        neutral = self.__m_blendshapes.neutral_pose()
        A = A[self.__m_marker_indices, :]
        w = np.linalg.solve( A.T@A, A.T@(marker_pose - neutral) )
        w = np.clip(w, a_min=0.0, a_max=1.0)
        return self.__m_blendshapes.make_pose_by_weight(w)
        

    def solve_phi_yt(self, Asums, bsums, B, prev_x_t_1, prev_x_acc):
        I = sp.identity(len(self.__m_sp_mass_matrix))
        h = self.__m_step_size
        h2 = self.__m_step_size**2
        M_inv = (1.0/self.__m_sp_mass_matrix)
        tmp1 = inv(I -  h2 * M_inv@Asums)
        tmp2 = h2@M_inv @ B
        phi = tmp1 @ tmp2 
        

        yt = tmp1@(prev_x_t_1, h*prev_x_acc + h2*M_inv@bsums)

        return phi, yt

    
    def solve_ut(self, phi, yt, dt):
        S = self.__M_S_sp_mat

        S_Phi = S @ phi
        S_Phi_T_S_Phi = S_Phi.T @ S_Phi
        factor= spchol(S_Phi_T_S_Phi)
        result_u = factor(S_Phi.T@ (dt - S@yt))
        return result_u
        
    def __system_forces(self, x):
        
        self.__m_precomputed_Ai_list = []



        func_list = [functools.partial(self.__displacement_constraints_A, k_p = self.__m_kp , neutral_v = self.__m_blendshapes.neutral_pose())]
        func_list = [functools.partial(self.__stretching_constraint_A, k_p = self.__m_kp , neutral_v = self.__m_blendshapes.neutral_pose())]
        func_list = [functools.partial(self.__bend_constraint_A, k_p = self.__m_kp  sp_laplacian =  self.__m_sp_laplacian_matrix, neutral_vv =self.__m_blendshapes.neutral_pose() )]
        self.__stretching_constraint_A(data)
        
        self.__displacement_constraints_A(data, rows, cols, self.__m_kp, self.__m_blendshapes.neutral_pose())
        for func in func_list:
            data = []
            rows = []
            cols = []
            func(data, rows, cols)
            self.__m_precomputed_Ai_list.append(sp.csc_matrix((data, (rows, cols)), dtype=np.float64))



    


    def simulate_time_step(self, x_prev, x_acc_prev, u_t):

        expr_poses = self.__m_blendshapes.expression_pose()
        x_acc_t = x_acc_prev + self.__m_step_size*(1.0/self.__M_S_sp_mat)@(expr_poses@u_t + self.__system_forces(x_t))
        x_t = x_acc_t + x_prev 

        return x_t, x_acc_prev

    def update(self, new_marker_pos):
        if self.__m_first_iter_flag:
            self.__m_first_iter_flag = False 
            self.__x_prev = self.__static_solve(new_marker_pos)
            self.__x_acc_prev = 0 

        for _ in range(self.__m_iteration_num):
            x_t = self.__x_prev + self.__m_step_size* self.__x_acc_prev
            A_sums, b_sums = self.__linearlize_forces(x_t)

            if self.__m_As_sum is None :
                self.__m_As_sum = A_sums
                self.__m_bs_sum = b_sums
            else:
                self.__m_As_sum = A_sums
                self.__m_bs_sum = b_sums
            
            phi, y_t = self.solve_phi_yt(self.__m_As_sum, self.__m_bs_sum )
            u_t = self.solve_ut(phi, y_t, new_marker_pos)
            x_t, x_acc_t = self.simulate_time_step(self.__x_prev, self.__x_acc_prev, u_t)

        self.__x_prev, self.__x_acc_prev = x_t, x_acc_t
        return x_t


if __name__ == "__main__":
    import os , glob 
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
        


    bshapes = blendshapes.Blendshapes(neutral, bs_list)
    bb = BlendForces()
    bb.set_blendshapes(bld=bshapes)
    bb.add_marker_index()
    
    bb.precompute()



    
    xt = bb.update()

    