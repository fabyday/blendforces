import numpy as np 





orig = np.random.rand(3,3)
v = np.random.rand(3,3)
f = np.array([[0,1,2]])


import igl

# orig, f = igl.read_triangle_mesh("./data/tri.obj")
# v = orig+np.random.normal(size=(orig.shape))
def edge_distance(v, f):
    res = [] 
    for f_list in f:
        v1 = f_list[0]
        v2 = f_list[1]
        v3 = f_list[2]

        e1 = v[v1] - v[v2]
        
        e2 = v[v2] - v[v3]
        e3 = v[v3] - v[v1]
        res+= [np.linalg.norm(e) for e in [e1,e2,e3]]
    return np.array(res)


wi = 1
A = np.zeros((len(v) * 3,len(v) * 3), dtype=np.float64)

b = np.zeros((len(v)*3, 1))



rest_pose_dist = edge_distance(orig, f)


xi = 0
yi = 1
zi = 2


def add_index(A, fr, ft):
    A[fr*3 + 0, fr*3 + 0 ] += 1
    A[fr*3 + 1, fr*3 + 1 ] += 1
    A[fr*3 + 2, fr*3 + 2 ] += 1
    A[fr*3 + 0, ft*3 + 0 ] += -1
    A[fr*3 + 1, ft*3 + 1 ] += -1
    A[fr*3 + 2, ft*3 + 2 ] += -1
    
    A[ft*3 + 0, fr*3 + 0 ] += 1
    A[ft*3 + 1, fr*3 + 1 ] += 1
    A[ft*3 + 2, fr*3 + 2 ] += 1
    A[ft*3 + 0, fr*3 + 0 ] += -1
    A[ft*3 + 1, fr*3 + 1 ] += -1
    A[ft*3 + 2, fr*3 + 2 ] += -1


for fi, f_list in enumerate(f):
    f_offset = 9*fi
    edge_offset = 3
    v1=  f_list[0]
    v2 =  f_list[1]
    v3 = f_list[2]


    e0 = np.array([v1,v2])
    e2 =  np.array([v2,v3])
    e3 =  np.array([v3,v1])

    fr, ft = e0
    add_index(A, fr, ft)
    add_index(A, ft, fr)
    
    fr, ft = e2
    add_index(A, fr, ft)
    add_index(A, ft, fr)

    fr, ft = e3
    add_index(A, fr, ft)
    add_index(A, ft, fr)
    
    distances = edge_distance(v,f_list.reshape(1,3))
    rest_pose_dist = edge_distance(orig,f_list.reshape(1,3))
    norm_r = rest_pose_dist.reshape(-1,1)/distances.reshape(-1,1)


    # # A = A.T + A
    # # A[list(range(len(A))), list(range(len(A)))] = 1
    # s = A.T[np.triu_indices(len(A), +1)]
    # A[np.tril_indices(len(A), -1)] += A[np.triu_indices(len(A), 1)]



    v1_off = 3*v1
    v2_off = 3*v2
    v3_off = 3*v3
    b[v1_off:v1_off+3, : ] += norm_r[0] * (v[v1] - v[v2]).reshape(3,1)
    b[v2_off:v2_off+3, : ] += norm_r[0] * (v[v2] - v[v1]).reshape(3,1)
    b[v2_off:v2_off+3, : ] += norm_r[1] * (v[v2] - v[v3]).reshape(3,1)
    b[v3_off:v3_off+3, : ] += norm_r[1] * (v[v3] - v[v2]).reshape(3,1)
    b[v3_off:v3_off+3, : ] += norm_r[2] * (v[v3] - v[v1]).reshape(3,1)
    b[v1_off:v1_off+3, : ] += norm_r[2] * (v[v1] - v[v3]).reshape(3,1)


np.set_printoptions(precision=2,suppress=True)

atb = A.T@b
print(atb)
a, s = np.linalg.eigh(A)
print(a)
print(s)
A = A + np.identity(len(A))
Ata = A.T@A 
aTAinv = np.linalg.inv(A.T@A)
bb = b + orig.reshape(-1,1)
print((aTAinv@A.T@bb).reshape(-1,3))

x = np.linalg.lstsq(A, b)
print(x[0].reshape(-1,3))
print(v)
print(orig)
print("det", np.linalg.det(A))
x = np.linalg.solve(Ata, atb)      

import matplotlib as plt 




print(x)

