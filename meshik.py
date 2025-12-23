import mesh 

import geometry_helper as geo 
import scipy.sparse as sp
from scipy.linalg import polar
import scipy 
# //1
# xi = asdasd
# for 1-10
# 	w = M(w+delta)
# 	A,b = linearize()
# 	x, xt = timeStep(w,A,b)
# //2
# for 1-10 --->meshik forloop
# 	w = M(w+delta)
# for 1-10 ---> for loop
# (Gx - T(w))
import numpy as np 
class FeatureVector:
    def __init__(self):
        pass


class MeshIKConstraint:
    def __init__(self):
        
        pass 
    
    
    def add_ref(self, ref : mesh.Mesh ):
        self.m_ref = ref
        
        
    def add_examples(self, examples ):
        self.m_examples = examples

    
    def _compute_extended_v(self, v):
        vn = geo.compute_face_normals(v, self.m_ref.f)
        vn += v[self.m_ref.f[:, 0], :]
        extended_v = np.concatenate([v, vn], 0)
        return extended_v, vn
    
    
    def __precompute_linearized_tw_B(self):
        d = self.mat_G@self.extended_v.reshape(-1,1)
        res = np.zeros((d.shape[0], len(self.m_examples_extended_v)))        
        for i in range(len(self.m_examples_extended_v)):
            res[:,i] = ( (self.mat_G @ self.m_examples_extended_v[i].reshape(-1,1)) - d ).reshape(-1)
        return res
    
    def get_reduced_linearized_tw_B(self, idx):
        G, c, _ = self.getCMatrixAndReducedG(idx)
        d = G @ self.extended_v.reshape(-1,1)
        res = np.zeros( (d.shape[0], len(self.m_examples_extended_v))        )
        for i in range(len(self.m_examples_extended_v)):
            res[:,i] = ( (G @ self.m_examples_extended_v[i].reshape(-1, 1)) - d ).reshape(-1)
        return res
        
    def precompute(self, compute_weigthed = 1.0):
        print("precompute")
        if not hasattr(self, "m_ref"):
            raise  ValueError("Reference mesh is not set.")
        self.extended_v, self.vn = self._compute_extended_v(self.m_ref.v)
        new_range = len(self.m_ref.v)+np.arange(len(self.vn)).reshape(-1,1)
        self.extended_f = np.concatenate([self.m_ref.f,  new_range],-1)
        
        
        self.m_examples_extended_v = []
        for i in range(len(self.m_examples)):
            self.m_examples_extended_v.append(self._compute_extended_v(self.m_examples[i].v)[0])
        
        
        self.ref_T = self._compute_T()
        self.inv_ref_T = self._compute_invT(self.ref_T)
        self.mat_G = self._compute_G(self.inv_ref_T)

        
        self.examples_T_list = self._computeExamplesT()
        self.Us = []
        self.Ps = []
        self.logUs = []
        for i in range(len(self.examples_T_list)):
            log_R_tmp = np.zeros((len(self.examples_T_list[i]), 3))
            U, P = self.RS_decompose(self.examples_T_list[i])
            for i in range(len(U)):
                log_R_tmp[i, ...] = self._rodrigues_log(U[i])
            self.Us.append(U)
            self.Ps.append(P)
            self.logUs.append(log_R_tmp)
        
        self.__B = self.__precompute_linearized_tw_B()
        print("precompute done")
    
    
    @property
    def B(self)    :
        return self.__B
    
    def _computeExamplesT(self):
        examples_T_list = []
        for ext_v in self.m_examples_extended_v:
            T = self.mat_G @ ext_v.reshape(-1,1)
            T=T.reshape(-1,3,3)
            examples_T_list.append(T)
        return examples_T_list        
    
        
        
    def _compute_G(self, ref_T):
        rows = 9*len(self.extended_f)
        cols = 3*len(self.extended_v)
        row_data = []
        col_data = []
        data = []
        
        for i, f in enumerate(self.extended_f):
            vi1,vi2,vi3, vi4 = f
            ref_T_i = ref_T[i]
            
            for dim in range(3):
                for k in range(3):
                    a = ref_T_i[:, k][0]
                    d = ref_T_i[:, k][1]
                    g = ref_T_i[:, k][2]
                    row_data.extend([i*9 + 3*dim + k] * 6)
                    col_data.extend([vi1*3+dim, vi2*3+dim, vi3*3+dim, vi4*3+dim, vi4*3+dim, vi4*3+dim])
                    data.extend([a, d,g, -a,-d,-g])
                    
        return sp.coo_matrix((data, (row_data, col_data)), shape=(rows, cols) ).tocsc()
    
    def _compute_invT(self, T):
        return np.linalg.inv(T)
    def _compute_T(self):
        
        v14 = self.extended_v[self.extended_f[:, 0], :] - self.extended_v[self.extended_f[:, 3], :]
        v24 = self.extended_v[self.extended_f[:, 1], :] - self.extended_v[self.extended_f[:, 3], :]
        v34 = self.extended_v[self.extended_f[:, 2], :] - self.extended_v[self.extended_f[:, 3], :]
        
        # F x (3x3) matrix
        res = np.hstack([v14,v24,v34])
        # convert rowwise vector to columnwise vector
        return np.transpose(res.reshape(-1,3,3), [0,2,1])
    
    
    def _exp_rotation(self, logR):
        """ compute the rotation matrix from the Rodrigues' rotation vector."""
        theta = np.linalg.norm(logR)
        
        
        
        if theta < 1e-8:
            return np.eye(3)
        
        axis = logR / theta
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        
        skew_symmetric = np.array([[0, -axis[2], axis[1]],
                                   [axis[2], 0, -axis[0]],
                                   [-axis[1], axis[0], 0]])
        
        R = cos_theta * np.eye(3) + sin_theta * skew_symmetric + (1 - cos_theta) * np.outer(axis, axis)
        return R
    
    def _rodrigues_log(self, R):
        """ compute the Rodrigues' rotation vector from a rotation matrix R."""
        eps = 1e-8
        trace = np.trace(R)
        cos_theta = (trace - 1) / 2.0
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        theta = np.arccos(cos_theta)

        if theta < eps:
            return np.zeros(3)

        skew = (R - R.T) / ( 2 * np.sin(theta))
        rx = skew[2,1]
        ry = skew[0,2]
        rz = skew[1,0]
        return theta * np.array([rx, ry, rz])
    
    def _roation_projection(self, ws, logRList):
        
        avg_logR = np.zeros(3, dtype=np.float64)
        for i in range(len(logRList)):
            avg_logR += logRList[i]*ws[i]
            
        return self._exp_rotation(avg_logR)
    
    
    def _scale_projection(self, ws, Si_list):
        avg_S = np.zeros((3,3), dtype=np.float64)
        for i in range(len(Si_list)):
            avg_S += Si_list[i] * ws[i]
        
        return avg_S
    
    def TwLinear(self, w):
        neww = w.reshape(-1)
        d = self.mat_G@self.extended_v.reshape(-1,1)
        # res = np.zeros_like(d)        
        # for i in range(len(self.m_examples_extended_v)):
        #     res += ( (self.mat_G @ self.m_examples_extended_v[i].reshape(-1,1)) - d )*neww[i]
        # return d + res.reshape(-1,1)
        return d+self.B@neww.reshape(-1,1)
    
    def Tw(self, w):
        size_t = len(self.logUs[0])
        result = np.zeros_like(self.Ps[0])
        identity = np.eye(3,3, dtype=np.float64)
        for i in range(size_t):
            R = self._roation_projection(w, [self.logUs[j][i] for j in range(len(self.logUs))])
            S = self._scale_projection(w, [self.Ps[j][i] for j in range(len(self.Ps))])
            
            result[i, :, :] = R @ (S + identity*(1 - w.sum()))
        return result
    def getCMatrixAndReducedG(self, idx):
        """
        return reducedG, c matrix
        """
        colidx = [] 
        for ii in idx:
            colidx.extend([ii*3,ii*3+1,ii*3+2])
        # colidx.sort()
        # COO로 변환
        coo = self.GMatrix().tocoo()
        # 제외할 컬럼만 0으로 만들기 위해 mask 사용
        mask = ~np.isin(coo.col, colidx)
        # 새로운 COO 생성
        new_coo = sp.coo_matrix((coo.data[mask], (coo.row[mask], coo.col[mask])),
                            shape=coo.shape)

        # 다시 CSC로 변환
        reduced_G  = new_coo.tocsc()
        Gmat = self.GMatrix()
        all_cols = np.arange(Gmat.shape[1])
        keep_cols = np.setdiff1d(all_cols, colidx)
        
        excluded_colsG = Gmat[:, keep_cols]
        return reduced_G, Gmat[:, colidx], excluded_colsG
    
            
    def RS_decompose(self, T):
        """
        Decompose the transformation matrix T into rotation and scaling.
        """
        Us = np.zeros_like(T)
        Ps = np.zeros_like(T)
        for i in range(len(T)):
            u, p = polar(T[i])
            
            if not np.allclose([u[:, 0, None].T@u[:,1], u[:, 0, np.newaxis].T@u[:,-1],  u[:, 0, None].T@u[:,1]], [0.0, 0.0,0.0]):
                print("Warning: Rotation matrix is not orthogonal.")
                
            Us[i,:,:]= u
            Ps[i:,:] = p
        
        return Us, Ps # Us: Rotations, Ps: scaling matrixs
    
    def gradientTw(self, w):
        wsize = w.size
        size_t = len(self.logUs[0])
        result = np.zeros_like(self.Ps)
        identity = np.eye(3,3, dtype=np.float64)
        
        def vec2Skew(r):
            rx, ry, rz = r
            return np.array([
                [0,   -rz,  ry],
                [rz,   0,  -rx],
                [-ry, rx,   0]
            ])
        
        for k in range(wsize):
            for i in range(size_t):
                R = self._roation_projection(w, [self.logUs[j][i] for j in range(len(self.logUs))])
                Rkj = self.logUs[k][i]
                
                S = self._scale_projection(w, [self.Ps[j][i] for j in range(len(self.Ps))])
                
                result[k,i, :, :] = R @ vec2Skew(Rkj) @ (S + identity*(1 - w.sum())) + R @ self.Ps[k][i]
        result = result.reshape(wsize, -1)
        return result
        
        
    
    
    
    def compute_G_from_two_vector(self):
        """this make G matrix with two vector per faces. normal vector won't be included."""
        T = self._compute_T_from_two_vector( )
        invT = self._compute_invT_from_two_vector(T)
        return self._compute_G_from_two_vector(invT)
        
    def _compute_invT_from_two_vector(self, T):
        invT = np.zeros_like(T)
        invT = np.transpose(invT, [0,2,1]) 
        # 3x2
        for i in range(len(T)):
            q,r = np.linalg.qr(T[i])
            #3 x 3, 3 x 2 
            Rj = r[:2, :2]
            Qja = q[:, :2]
            invT[i, :, :] = np.linalg.inv(Rj)@Qja.T
        return invT
    
    def _compute_T_from_two_vector(self):
        
        v14 = self.m_ref.v[self.m_ref.f[:, 1], :] - self.m_ref.v[self.m_ref.f[:, 0], :]
        v24 = self.m_ref.v[self.m_ref.f[:, 2], :] - self.m_ref.v[self.m_ref.f[:, 0], :]
        
        # F x (3x3) matrix
        res = np.hstack([v14,v24])
        return np.transpose(res.reshape(-1,2,3), [0,2,1]) # N,3,2
        
        
        # convert rowwise vector to columnwise vector
        return
    
    def _compute_G_from_two_vector(self, ref_T):
        rows = 9*len(self.m_ref.f)
        cols = 3*len(self.m_ref.v)
        row_data = []
        col_data = []
        data = []
        
        for i, f in enumerate(self.m_ref.f):
            vi1,vi2,vi3 = f
            ref_T_i = ref_T[i]
            
            for dim in range(3):
                for k in range(3):
                    a = ref_T_i[:, k][0]
                    d = ref_T_i[:, k][1]
                    row_data.extend([i*9 + 3*dim + k] * 4)
                    col_data.extend([vi1*3+dim, vi2*3+dim, vi3*3+dim, vi3*3+dim])
                    data.extend([a,d,-a,-d])
                    
        return sp.coo_matrix((data, (row_data, col_data)), shape=(rows, cols) ).tocsc()
        
        
    def GMatrix(self):
        return self.mat_G
    
    
    def GMatrix_with_v4(self):
        return self.mat_G[:, :3*len(self.m_ref.v)]
    
    
    
    def solve_nonlinear_from_lmk(self, reudcedG, cMatrix, lmk, lmk_indces, init_w = None, init_x = None, iter_num = 10):
        c = cMatrix @ lmk.reshape(-1,1)
        w = None 
        if init_w is None:
            w = np.zeros((len(self.m_examples_extended_v), 1))
        else:
            w = init_w.reshape(-1,1)        
        oA = reudcedG 
        # A = A + 1e-8*sp.identity(A.shape[0])
        for i in range(iter_num):
            print("  iter ", i)
            JJJ = self.gradientTw(w)
            J = sp.csc_array(JJJ)
            A = sp.hstack([oA, -J.T])
            AtA = A.T@A + 1e-8*sp.identity(A.shape[-1])
            Mw = self.Tw(w).reshape(-1,1)
            
            b = A.T@(Mw-c)
            xt = sp.linalg.spsolve(AtA ,b)
            dsize = len(w)
            delta = xt[-dsize:, None]
            w = w + delta
            
        ret_xt = np.zeros_like(self.m_ref.v)
        mask = np.setdiff1d(np.arange(len(self.m_ref.v)), lmk_indces)
        xt = xt[:-dsize].reshape(-1,3)
        nn_xt =xt[:-len(self.m_ref.f), :]
        ret_xt[mask ,  : ] = nn_xt
        ret_xt[lmk_indces, :] = lmk
        return ret_xt
    
    def solve(self,Tw):
        return sp.linalg.spsolve(self.mat_G.T@self.mat_G, self.mat_G.T@Tw.reshape(-1,1)).reshape(-1,3)
    
    
    def solve_from_lmk_linear(self, reudcedG, cMatrix, lmk, lmk_indces):
        c = cMatrix @ lmk.reshape(-1,1)
        w = np.zeros((len(self.m_examples_extended_v), 1))
        A = reudcedG 
        # A = A + 1e-8*sp.identity(A.shape[0])
        AtA = A.T@A + 1e-8*sp.identity(A.shape[-1])
        b = A.T@(self.TwLinear(w)-c)
        xt = sp.linalg.spsolve(AtA ,b)
        d = self.mat_G@self.extended_v.reshape(-1,1)
        w = sp.linalg.spsolve(self.B.T@self.B, self.B.T@((A @xt.reshape(-1,1) + c) - d.reshape(-1,1)) )
        xt = sp.linalg.spsolve(AtA, A.T@(self.TwLinear(w)-c)).reshape(-1,3)
        xt[lmk_indces, :] = lmk
        return xt
        
        
    
    
if __name__ == "__main__":
    # Example usage
    import os, glob
    import mesh as mm
    meshik = MeshIKConstraint()
    
    
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
        
    lmk_idx = [1278,1272,12,1834,243,781,2199,1447,966,3661,4390,3022,2484,4036,2253,3490,3496,268,493,1914,2044,1401,3615,4240,4114,2734,2509,978,4527,4942,4857,1140,2075,1147,4269,3360,1507,1542,1537,1528,1518,1511,3742,3751,3756,3721,3725,3732,5708,5695,2081,0,4275,6200,6213,6346,6461,5518,5957,5841,5702,5711,5533,6216,6207,6470,5517,5966,]

    meshik.add_ref(neutral)
    meshik.add_examples(bs_list)
    meshik.precompute()
    
    
    m_list = []
    
    adder = 0.1
    # adder = 0.5
    ws = np.zeros(len(bs_list))
    for i in range(len(bs_list)):
        while ws[i] < 1.0:
            ws[i] += 2*adder
            ws -= adder
            ws = np.clip(ws, 0.0, 1.0)
            result = meshik.Tw(ws)
            # result = meshik.TwLinear(ws)
            print(f"{i} th / {len(bs_list)} and ws: {ws[i]}/1.0")
        
            m = meshik.solve(result)
            m = m[:len(meshik.m_ref.v), :]
            m_list.append(m)
            
            
    np.save("testanim_out/testik.npy", np.array(m_list))
    m_list = np.load("testanim_out/testik.npy")
    import viewer ,sys
    app = viewer.QApplication(sys.argv)

    win = viewer.MainWindow(neutral_pth, pipe=False, framerate=60)
    win.add_animation(m_list)
    
    win.show()
    sys.exit(app.exec_())
    print("ex")