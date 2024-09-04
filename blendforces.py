import numpy as np 
from . import mesh as mm
import typing 
import blendshapes
import scipy.sparse as  sp 
import geometry_helper as geo
from scipy.sparse.linalg import inv

stiffnesses_args = tuple | typing.Literal["auto"] 
mass_args = typing.Union[typing.List[float], int, np.ndarray]
class BlendForces:
    def __init__(self, mesh : mm.Mesh, stiffnesses  : stiffnesses_args = "auto", iteration_num : int = 10, step_size = 0.16, mass : mass_args = 10.0):
        self.__m_stiffnesses = stiffnesses
        self.__m_mesh = mesh 
        self.set_iteration_num(iteration_num)
        self.set_step_size(step_size)
        self.reset_simulation()
        self.__m_mass = mass

        self.__m_blendshapes : blendshapes.Blendshapes = None 


        #TODO for testing
        self.__m_kp = 1.0
        self.__m_ks = 1.0
        self.__m_kb = 1.0


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




    def __displacement_constraints(self, data,  row_indices, col_indices, k_p, neutral_v, v):
        assert(len(v) == len(neutral_v), "len() size between neutral v and v is diff")
        new_ind_row = row_indices[-1] + 1
        data += [1.0 for _ in range(len(neutral_v))]
        row_indices += [new_ind_row + i for i in range(len(neutral_v))]
        col_indices += [ i for i in range(len(neutral_v))]

        return data, row_indices, col_indices

            
        

    def __bend_constraint(self, data : list, row_indices : list, col_indices : list, sp_laplacian : sp.csc_matrix, neutral_vv ):
        row, col = sp_laplacian.nonzero()
        dt = sp_laplacian.data


        row_offset =row_indices[-1] + 1
        row += row_offset 
        data.extend(dt)
        row_indices.extend(row)
        col_indices.extend(col)
        return 
        

    def __stretching_constraint(self, data :list, row_indices : list, col_indices : list, v : np.ndarray, edges, neutral_rest_length ):
        new_row_offset = row_indices[-1] + 1
        data += len(edges)*[1,-1]

        for i, j in edges:
            row_indices.append(new_row_offset)
            col_indices.append(i)
            col_indices.append(j)
            new_row_offset += 1


    







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



    def add_marker_index(self, index_list : list[int]):
        self.__m_marker_indices = index_list
        data = [1 for _ in len(self.__m_marker_indices)]
        row = [ i for i in range(len(self.__m_marker_indices))]
        col = self.__m_marker_indices

        N = len(self.__m_blendshapes.neutral_pose())
        self.__M_S_sp_mat = sp.csc_matrix((data, (row,col)), shape=(len(self.__m_marker_indices),N), dtype=np.float64) 


    def __linearlize_forces(self, x_t):
        data = [] 
        row_indices = []
        col_indices = []

        self.__displacement_constraints(data,  row_indices, col_indices, self.__m_blendshapes.neutral_pose(), x_t)
        self.__stretching_constraint(data,  row_indices, col_indices)
        self.__bend_constraint(data,  row_indices, col_indices)
        


        sp.csc_matrix((data, (row_indices, col_indices)), shape=())


    def static_solve(self, marker_pose):
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
        pseudo_S_phi_inv = inv(S_Phi.T @ S_Phi)
        return pseudo_S_phi_inv @ (dt - S@yt)
        
    def __system_forces(self, x):
        pass
    


    def simulate_time_step(self, x_prev, x_acc_prev, u_t):

        expr_poses = self.__m_blendshapes.expression_pose()
        x_acc_t = x_acc_prev + self.__m_step_size*(1.0/self.__M_S_sp_mat)@(expr_poses@u_t + self.__system_forces(x_t))
        x_t = x_acc_t + x_prev 

        return x_t, x_acc_prev

    def update(self, new_marker_pos):
        if self.__m_first_iter_flag:
            self.__m_first_iter_flag = False 
            self.__x_prev = self.static_solve(new_marker_pos)
            self.__x_acc_prev = 0 


        for _ in range(self.__m_iteration_num):
            x_t = self.__x_prev + self.__m_step_size* self.__x_acc_prev
            A, b = self.__linearlize_forces(x_t)

            if self.__m_As_sum is None :
                self.__m_As_sum += A
                self.__m_bs_sum += b
            else:
                self.__m_As_sum = A
                self.__m_bs_sum = b 
            
            phi, y_t = self.solve_phi_yt(self.__m_As_sum, self.__m_bs_sum )
            u_t = self.solve_ut(phi, y_t, new_marker_pos)
            x_t, x_acc_t = self.simulate_time_step(self.__x_prev, self.__x_acc_prev, u_t)

        self.__x_prev, self.__x_acc_prev = x_t, x_acc_t
        return x_t





