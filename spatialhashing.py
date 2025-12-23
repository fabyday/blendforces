import numpy as np
from collections import defaultdict
import mesh as mm
import geometry_helper as geo 
import scipy.sparse as sp

from numpy.polynomial.polynomial import Polynomial


class OptimSpatialHashGrid:
    # https://matthias-research.github.io/pages/publications/tetraederCollision.pdf
    
    def __init__(self):
        self.tau = 0.001
        self.cell_size = 0.1
        self.hashtable_size = 10000
        self.hash_table = defaultdict(lambda : defaultdict)
        self.p1 = 73856093
        self.p2 = 19349663
        self.p3 = 83492791
        self.current_timestamp = 0

    def _hash(self, point):
        point = np.floor(point / self.cell_size).astype(int)
        return (point[0]*self.p1 ^ point[1]*self.p2 ^ point[2]*self.p3)% self.hashtable_size
    
    def ee_collision(self, x0, x1, x2, x3):
        w0 = 1
        w1 = 1
        w2 = 1 
        w3 = 1
        return x0*w0 + w1*x1 - w2*x2 - w3*x3
        
    def vt_collision(self, x0, x1, x2, x3):
        w0 = 1
        w1 = 1
        w2 = 1 
        w3 = 1
        return x0*w0 + w1*x1 + w2*x2 - w3*x3
    
    def _compute_average_edge_length(self , V, F):
        # 모든 엣지를 추출 (삼각형의 세 변)
        edges = np.concatenate([
            F[:, [0, 1]],
            F[:, [1, 2]],
            F[:, [2, 0]]
        ], axis=0)

        # 정렬하여 중복 제거 (예: [2, 5]와 [5, 2]를 동일하게 취급)
        edges = np.sort(edges, axis=1)
        edges = np.unique(edges, axis=0)

        # 각 엣지의 길이 계산
        edge_vectors = V[edges[:, 0]] - V[edges[:, 1]]
        edge_lengths = np.linalg.norm(edge_vectors, axis=1)

        # 평균 길이 반환
        return np.mean(edge_lengths)
    
    # def firststep_test(self, triid1, triid2, data, rows, cols):
    #     # The first test computes barycentric coordinates of 
    #     # a vertex with respect to a vertex of a tetrahedron
    #     t1_vindices = self.f[triid1]
    #     t2_vindices = self.f[triid2]
    #     v_tri_candidates = [] 
    #     for t2vi in t2_vindices:
    #         if t2vi in t1_vindices:
    #             continue 
    #         v_tri_candidates.append(t2vi)
    #     result = []
    #     v1,v2,v3 = t1_vindices
    #     e1 = self.ref_v[v1] - self.ref_v[v3]
    #     e2 = self.ref_v[v2] - self.ref_v[v3]
    #     normal = np.cross(e1,e2)
    #     if np.linalg.norm(normal) < 1e-8:
    #         return  # degenerate triangle
        
        
    #     normalized_normal =  normal/np.linalg.norm(normal)
    #     e1e1d = np.dot(e1,e1)
    #     e1e2d = np.dot(e1,e2)
    #     e2e2d = np.dot(e2, e2)
    #     A = np.array([[e1e1d, e1e2d],[e1e2d, e2e2d]])
    #     Ainv = np.linalg.inv(A)
    #     for v2i in v_tri_candidates:
    #         v43 = self.ref_v[v2i] - self.ref_v[v3]
            
    #         s = np.abs(v43@normalized_normal)
    #         if s < 0.01:
                
                
            
    #             b= np.array([[np.dot(e1, v43)],
    #             [ np.dot(e2, v43)]])
    #             w = Ainv@b
    #             w1, w2, w3 = w[0,0], w[1,0], 1 - (w[0,0] + w[1,0])
    #             if w1 < 0 or w2 < 0  or w1 > 1 or w2 > 1 or w3<0 or w3 > 1:
    #                 continue
                
    #             ns = (w1*self.ref_v[v1] + w2*self.ref_v[v2] + w3*self.ref_v[v3] - 1*self.ref_v[v2i])@normal
                
    #             if ns > 0: # tau == 0.01mm
    #                 continue 
                
                    
    #             w1n = w1*normal 
    #             w2n = w2*normal 
    #             w3n = w3*normal
    #             w4n = -1*normal
    #             row_idx = len(data)//12
    #             data.append(w1n[0]);data.append(w1n[1]);data.append(w1n[2])
    #             data.append(w2n[0]);data.append(w2n[1]);data.append(w2n[2])
    #             data.append(w3n[0]);data.append(w3n[1]);data.append(w3n[2])
    #             data.append(w4n[0]);data.append(w4n[1]);data.append(w4n[2])
    #             rows.append(row_idx);rows.append(row_idx);rows.append(row_idx)
    #             rows.append(row_idx);rows.append(row_idx);rows.append(row_idx)
    #             rows.append(row_idx);rows.append(row_idx);rows.append(row_idx)
    #             rows.append(row_idx);rows.append(row_idx);rows.append(row_idx)
    #             cols.append(v1*3 + 0);cols.append(v1*3 + 1);cols.append(v1*3 + 2)
    #             cols.append(v2*3 + 0);cols.append(v2*3 + 1);cols.append(v2*3 + 2)
    #             cols.append(v3*3 + 0);cols.append(v3*3 + 1);cols.append(v3*3 + 2)
    #             cols.append(v2i*3 + 0);cols.append(v2i*3 + 1);cols.append(v2i*3 + 2)
                
    #             result.append([[v1, v2, v3, v2i], [w1,w2, w3], normal])
                
    #     return result
            
    def firststep_test(self, triid1, triid2, data, rows, cols):
        t1_vindices = self.f[triid1]
        
        t2_vindices = self.f[triid2]
        v_tri_candidates = []
        for t2vi in t2_vindices:
            if t2vi in t1_vindices:
                continue
            v_tri_candidates.append(t2vi)

        result = []
        v1, v2, v3 = t1_vindices
        e1 = self.ref_v[v1] - self.ref_v[v3]
        e2 = self.ref_v[v2] - self.ref_v[v3]
        normal = np.cross(e1, e2)
        norm_len = np.linalg.norm(normal)
        if norm_len < 1e-8:
            return  # degenerate triangle

        normal = normal / (norm_len +1e-8)

        # 삼각형 중심 계산 (normal 방향 보정용)
        tri_center = (self.ref_v[v1] + self.ref_v[v2] + self.ref_v[v3]) / 3

        # 바리센트릭 계산 준비
        e1e1d = np.dot(e1, e1)
        e1e2d = np.dot(e1, e2)
        e2e2d = np.dot(e2, e2)
        A = np.array([[e1e1d, e1e2d], [e1e2d, e2e2d]])
        Ainv = np.linalg.inv(A)

        for v2i in v_tri_candidates:
            v43 = self.ref_v[v2i] - self.ref_v[v3]

            s = np.abs(np.dot(v43, normal))
            if s < 0.01:  # 가까운 평면 위 점만 처리
                # normal 방향 보정
                # to_test = self.ref_v[v2i] - tri_center
                # if np.dot(to_test, normal) > 0:
                    # normal = -normal

                b = np.array([[np.dot(e1, v43)], [np.dot(e2, v43)]])
                w = Ainv @ b
                w1, w2 = w[0, 0], w[1, 0]
                w3 = 1 - (w1 + w2)

                if not (0 <= w1 <= 1 and 0 <= w2 <= 1 and 0 <= w3 <= 1):
                # if not abs(1-(w1+w2+w3)) < 1e-8 and not (0 <= w1 <= 1 and 0 <= w2 <= 1 and 0 <= w3 <= 1):
                    continue
                
                # if abs(1-(w1+w2+w3)) > 1e-8:
                #     continue
                
                
                # ns = (w1 * self.ref_v[v1] + w2 * self.ref_v[v2] + w3 * self.ref_v[v3] - self.ref_v[v2i]) @ normal
                ns = (self.ref_v[v2i] - ( w1 * self.ref_v[v1] + w2 * self.ref_v[v2] + w3 * self.ref_v[v3]) ) @ normal
                if ns > self.tau:
                    continue  # 침투가 아님
                
                depth = max(0, -ns)
                tau0 = self.tau 
                alpha = 0.5
                adaptive_tau = tau0 + alpha * depth
                adaptive_tau = self.tau 
                # 힘 벡터 구성
                # vn = self.ref_normal_v[v2i]/ ( np.linalg.norm(self.ref_normal_v[v2i]) + 1e-8)
                w1n = -w1 * normal
                w2n = -w2 * normal
                w3n = -w3 * normal
                # w4n = -1 * normal
                w4n = 1 * normal
                row_idx = len(data) // 12
                
                data.extend([w1n[0], w1n[1], w1n[2],
                            w2n[0], w2n[1], w2n[2],
                            w3n[0], w3n[1], w3n[2],
                            w4n[0], w4n[1], w4n[2]])
                rows.extend([row_idx] * 12)
                cols.extend([
                    v1 * 3 + 0, v1 * 3 + 1, v1 * 3 + 2,
                    v2 * 3 + 0, v2 * 3 + 1, v2 * 3 + 2,
                    v3 * 3 + 0, v3 * 3 + 1, v3 * 3 + 2,
                    v2i * 3 + 0, v2i * 3 + 1, v2i * 3 + 2
                ])

                result.append([[v1, v2, v3, v2i], [w1, w2, w3], normal, adaptive_tau])
        return result

    def secondstep_test(self, tri_vid1, tri_vid2, data, rows, cols):
        # The second test considers oriented faces of a tetrahedron and checks,
        # whether a vertex is in the positive or negative half–space of these faces.
        tri1_edges = self.get_edges(self.f[tri_vid1])
        tri2_edges = self.get_edges(self.f[tri_vid2])
        
        time_delta = 1/24 # 24 fps
        
        edge_candidates = []
        for t1ei in tri1_edges:
            for t2ej in tri2_edges: # test same edges
                if (t1ei[0] == t2ej[0] and t1ei[1] == t2ej[1]) \
                or (t1ei[1] == t2ej[0] and t1ei[0] == t2ej[1]):
                    continue
                #test share vertex
                if t1ei[0] == t2ej[0] or t1ei[0] == t2ej[1] or\
                    t1ei[1] == t2ej[0] or t1ei[1] == t2ej[1] :
                        continue
                edge_candidates.append([t1ei, t2ej])
        result= []
        for e1, e2 in edge_candidates:
            v1= self.ref_v[e1[0]] 
            v2 = self.ref_v[e1[1]]
            v3= self.ref_v[e2[0]] 
            v4= self.ref_v[e2[1]] 
            x21 = v2 - v1 
            x43 = v4 - v3
            x31 = v3 - v1
            x41 = v4 - v1
            vel1 = self.velocity[e1[0]]
            vel2 = self.velocity[e1[1]]
            vel3 = self.velocity[e2[0]]
            vel4 = self.velocity[e2[1]]
            vel21 = vel2 - vel1
            vel41 = vel4 - vel1
            vel31 = vel3 - vel1
            vel43 = vel4 - vel3
            
            line_normal = np.cross(x21, x43)
            cross_val = np.linalg.norm(line_normal)

            # degenerate case, parelles
            if cross_val < 10e-6:
                v43 = v4-v3
                # v43n = np.linalg.norm(v43)
                # distance = np.linalg.norm(np.cross(v1-v3, v43n))
                # degenerate case, parelles
                continue
            
            A = np.array([[np.dot(x21, x21), -np.dot(x21, x43)],\
                [-np.dot(x21, x43), np.dot(x43, x43)]])
            bb  = np.array([[np.dot(x21, x31)],[-np.dot(x43,x31)]])
            
            ab = np.linalg.inv(A) @ bb
            ab = ab.reshape(-1)
            a,b = ab[0],ab[1]
            # degenerate case, closest points lie outside the segments
            if a < 0 or a > 1 or b < 0 or b > 1:
                continue  
            p1 = v1 + a * x21
            p2 = v3 + b * x43
            distance = np.linalg.norm(p1 - p2)
                
            # if distance > self.tau:
            #     continue
            
            # (1-a)*v1 + a*v2
            # (1-b)*v3 + b*v4
            
            coeffs = [0.0, 0.0, 0.0, 0.0]
            mat = np.zeros((3, 3, 3))
            mat[0, 1, 2] = mat[1, 2, 0] = mat[2, 0, 1] = 1
            mat[2, 1, 0] = mat[1, 0, 2] = mat[0, 2, 1] = -1
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        e = mat[i, j, k]
                        if(e == 0):
                            continue
                        coeffs[3] += e * vel21[i] * vel31[j] * vel41[k]
                        coeffs[2] += e * (x21[i]*vel31[j]*vel41[k] + vel21[i]*x31[j]*vel41[k] + vel21[i]*vel31[j]*x41[k])
                        coeffs[1] += e * (x21[i]*x31[j]*vel41[k] + x21[i]*vel31[j]*x41[k] + vel21[i]*x31[j]*x41[k])
                        coeffs[0] += e * x21[i]*x31[j]*x41[k]
            roots_r = Polynomial(coeffs).roots() # a+bx+bx**2
            roots = [r.real for r in roots_r if 0<= np.isreal(r) <= time_delta]
            roots.sort()

            collision_flag = False 
            for r in roots:
                if np.isreal(r) and 0 <= r.real <= time_delta:
                    t = r 
                    p1 = v1 + a * (x21 + vel21 * t)
                    p2 = v3 + b * (x43 +vel43*t)
                    leng = np.linalg.norm(p1 - p2)
                    if leng < self.tau:
                        collision_flag = True 
                        continue 
            
            if not collision_flag:
                continue 
            # x21 + v

            depth = max(0.0,self.tau - distance)
            tau_cur = self.tau + 0.5*depth
            a1,a2,a3,a4 = 1-a, a, 1-b, b
            # a1*v1 + a2*v2
            # a3*v3 + a4*v4
            # line_normal = (1-b)*self.ref_normal_v[e2[0]] + b*self.ref_normal_v[e2[1]]
            line_normal /= (np.linalg.norm(line_normal) + 1e-8)
            row_idx = len(data)//12
            
            e11d = line_normal*a1
            e12d = line_normal*a2
            e21d = -line_normal*a3    
            e22d = -line_normal*a4
            data.append(e11d[0]);data.append(e11d[1]);data.append(e11d[2])
            data.append(e12d[0]);data.append(e12d[1]);data.append(e12d[2])
            data.append(e21d[0]);data.append(e21d[1]);data.append(e21d[2])
            data.append(e22d[0]);data.append(e22d[1]);data.append(e22d[2])
            rows.append(row_idx);rows.append(row_idx);rows.append(row_idx)
            rows.append(row_idx);rows.append(row_idx);rows.append(row_idx)
            rows.append(row_idx);rows.append(row_idx);rows.append(row_idx)
            rows.append(row_idx);rows.append(row_idx);rows.append(row_idx)
            cols.append(3*e1[0]+0);cols.append(3*e1[0]+1);cols.append(3*e1[0]+2)
            cols.append(3*e1[1]+0);cols.append(3*e1[1]+1);cols.append(3*e1[1]+2)
            cols.append(3*e2[0]+0);cols.append(3*e2[0]+1);cols.append(3*e2[0]+2)
            cols.append(3*e2[1]+0);cols.append(3*e2[1]+1);cols.append(3*e2[1]+2)
            
            result.append([[e1[0],e1[1], e2[0], e2[1]], [a1,a2,a3,a4], line_normal, tau_cur])
            # (1-a)*v1 + a*v2
            # (1-b)*v3+b*v4
            # v3 + x43*b
        return result 
            
    
    def attach_face_data(self, fdata)    :
        """
         Nx3 tri face data
        """
        self.f = fdata 
    def upate_vertex_data(self, v):
        if hasattr(self, "ref_v"):
            self.prev_ref_v = self.ref_v
        self.ref_v = v
        if hasattr(self, "velocity"):
            time = 1/24
            self.velocity =  (self.ref_v - self.prev_ref_v) / time
        else : 
            self.velocity = np.zeros_like(v)
        
        self.ref_normal_v = geo.compute_vertex_normals(self.ref_v, self.f)
        if(self.current_timestamp == 0 ):
            self.cell_size = self._compute_average_edge_length(self.ref_v, self.f)
        
        # for vi in range(len(self.ref_v)):
            # self.insert_vertex(vi) 
        self.current_timestamp += 1
        for tri_id in range(len(self.f)):
            self.insert_triangle(tri_id)
        
    def _cell(self, xyz):
        return np.floor(xyz / self.cell_size).astype(int)
    
    def _hash_from_cell(self, cell_xyz):
        return (
        (int(cell_xyz[0]) * self.p1) ^
        (int(cell_xyz[1]) * self.p2) ^
        (int(cell_xyz[2]) * self.p3)
    ) % self.hashtable_size

        
    
    def _get_covered_cells(self, aabb_min, aabb_max):
        
        min_cell = self._cell(aabb_min)
        max_cell = self._cell(aabb_max)
        for x in range(min_cell[0], max_cell[0] + 1):
            for y in range(min_cell[1], max_cell[1] + 1):
                for z in range(min_cell[2], max_cell[2] + 1):
                    yield self._hash_from_cell(np.array([x, y, z]))
        
    def insert_vertex(self, vid):
        hashval = self._hash(self.ref_v[vid])
        self.hash_table[hashval].append(vid)
        
    def insert_triangle(self, tri_ind):
        fi_v_indices = self.f[tri_ind]
        vertices = self.ref_v[fi_v_indices]
        
        
        aabb_min = np.min(vertices, axis=0) - self.tau
        aabb_max = np.max(vertices, axis=0)  + self.tau
        for cell in self._get_covered_cells(aabb_min, aabb_max):
            if (cell not in self.hash_table) or (self.hash_table[cell]["timestamp"] != self.current_timestamp):
                self.hash_table[cell] = {"timestamp" : self.current_timestamp, "tri_ind" : []}
            self.hash_table[cell]["tri_ind"].append(tri_ind)
            


    def get_cells(self):
        return list(self.hash_table.keys())

    
    def _aabb_overlap(self, min1, max1, min2, max2):
        return np.all(max1 >= min2) and np.all(max2 >= min1)


    def query_overlapped_tris(self):
        # result = [[] for f in range(len(self.f))]
        # for tri_id in range(len(self.f)):
        #     fi_v_indices = self.f[tri_id]
        #     vertices = self.ref_v[fi_v_indices]

        
        #     aabb_min = np.min(vertices - self.tau, axis=0)
        #     aabb_max = np.max(vertices + self.tau, axis=0)
        #     for cell in self._get_covered_cells(aabb_min, aabb_max):
        #         if self.hash_table[cell] is not None:
        #             if self.hash_table[cell]["timestamp"] == self.current_timestamp:
        #                 result[tri_id].extend(filter(lambda x: x != tri_id, self.hash_table[cell]["tri_ind"]))

        #     for i in range(len(result)):
        #         result[i] = list(set(result[i]))
        
        # return result
        checked_pairs = set()
        result_pairs = []
        for tri_id in range(len(self.f)):
            fi_v_indices = self.f[tri_id]
            vertices = self.ref_v[fi_v_indices]
            aabb_min = np.min(vertices , axis=0)- self.tau
            aabb_max = np.max(vertices , axis=0)+ self.tau

            for cell in self._get_covered_cells(aabb_min, aabb_max):
                if cell not in self.hash_table:
                    continue
                if self.hash_table[cell]["timestamp"] != self.current_timestamp:
                    continue

                for other_tri in self.hash_table[cell]["tri_ind"]:
                    if other_tri == tri_id:
                        continue
                    pair = tuple(sorted((tri_id, other_tri)))
                    if pair in checked_pairs:
                        continue

                    # ▶ 추가: 실제 AABB 교차 검사
                    other_vertices = self.ref_v[self.f[other_tri]]
                    other_min = np.min(other_vertices - self.tau, axis=0)
                    other_max = np.max(other_vertices + self.tau, axis=0)
                    if not self._aabb_overlap(aabb_min, aabb_max, other_min, other_max):
                        continue

                    checked_pairs.add(pair)
                    result_pairs.append(pair)

        return result_pairs
        
        for tri_id in range(len(self.f)):
            fi_v_indices = self.f[tri_id]
            vertices = self.ref_v[fi_v_indices]
            aabb_min = np.min(vertices - self.tau, axis=0)
            aabb_max = np.max(vertices + self.tau, axis=0)

            for cell in self._get_covered_cells(aabb_min, aabb_max):
                if self.hash_table[cell]["timestamp"] != self.current_timestamp:
                    continue
                for other_tri in self.hash_table[cell]["tri_ind"]:
                    if other_tri == tri_id:
                        continue
                    pair = tuple(sorted((tri_id, other_tri)))
                    if pair in checked_pairs:
                        continue
                    checked_pairs.add(pair)
                    result_pairs.append(pair)

        return result_pairs  # [(tri_id1, tri_id2), ...]
    
    def get_edges(self, tri_index):
        return [[tri_index[0], tri_index[1]], [tri_index[1], tri_index[2]], [tri_index[2], tri_index[0]]]
                        

    def calc_forces(self, candidates_tris_idx):
        eterm = {'vt': [], 'ee': []}
        rows= []
        cols = []
        datas = []
        taus = []
        for fid in candidates_tris_idx:
            vt = self.firststep_test(fid[0], fid[1], datas, rows, cols)
            if vt: 
                eterm['vt'].extend(vt)
                taus.extend(map(lambda x : x[-1], vt))
            
            ee = self.secondstep_test(fid[0], fid[1], datas, rows, cols)
            if ee: 
                eterm['ee'].extend(ee)
                taus.extend(map(lambda x : x[-1], ee))
        print(len(datas)//12, "vv:", len(eterm['vt']), "ee", len(eterm['ee']) )
        # 0.01*np.ones((len(datas)//12, 1))
        return eterm, sp.coo_matrix((datas, (rows, cols)), shape=( len(datas)//12 ,len(self.ref_v)*3)).tocsc(), np.array(taus).reshape(len(datas)//12, 1)

if __name__ == "__main__":

    import igl, os , glob
    import mesh as mm

    data_path = "D:\\lab\\2022\\mycode\\FaceCaptureWithIK\\data\\ICT-data"
    neutral_pth = os.path.join(data_path, "generic_neutral_mesh.obj")
    # neutral_pth = os.path.join("./data/test2.obj")
    # neutral_pth = os.path.join("./data/test3.obj")
    neutral_pth = os.path.join("./data/test5.obj")
    neutral = mm.Mesh()
    neutral.load_from_file(neutral_pth)

    opt = OptimSpatialHashGrid()
    opt.attach_face_data(neutral.f)
    import time
    import viewer
    start = time.time()
    opt.upate_vertex_data(neutral.v)
    end = time.time()
    print(f'{end - start} : vertex update time')

    start = time.time()
    a = opt.query_overlapped_tris()
    end = time.time()
    print(f'{end - start} : candidates query time')

    start = time.time()
    f = opt.calc_forces(a)
    end = time.time()
    print(f'{end - start} : force calctime')
    
    print("end")