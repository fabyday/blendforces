import igl 
from sksparse.cholmod import cholesky as spchol


def cot(pt):
    v, f = igl.read_triangle_mesh(pt)
    c = igl.cotmatrix(v,f)
    return c
    

v, f = igl.read_triangle_mesh("./data/quad.obj")
c = igl.cotmatrix(v,f)
print(c.toarray())
spchol(c)

s = spchol(cot("./data/test.obj"))
a = spchol(cot("./data/tri.obj"))
import os  
data_path = "D:\\lab\\2022\\mycode\\FaceCaptureWithIK\\data\\ICT-data"
neutral_pth = os.path.join(data_path, "generic_neutral_mesh.obj")


sss = -cot(neutral_pth)
sss
print(sss.getrow(0))
cols = sss.indices           # 0이 아닌 값이 있는 열의 인덱스
vals = sss.data              # 그 열의 실제 값
d = 0
for col, val in zip(cols[1:], vals[1:]):
    d += val
print(d)
print(vals[0])
print("wait")
v, f = igl.read_triangle_mesh(neutral_pth)
a = igl.doublearea(v,f)
import numpy as np 
print(np.any(a < 1e-12))

print(a.min())
diff = np.linalg.norm(sss - sss.T)
import scipy.sparse as sp
eps = 1e-6
sss += eps * sp.eye(sss.shape[0])  # εI 더해서 definite하게 만듬
spchol(sss)
print(diff)



coo = sss.tocoo()
row = coo.row 
col = coo.col 
d = coo.data
rows, cols = coo.shape
new_row = [] 
new_col = [] 
new_data = []
for r, c, dd in zip(row, col, d):
    for axis in range(3):
        new_row.append(3 * r + axis)
        new_col.append(3 * c + axis)
        new_data.append(dd)
        
re = sp.coo_matrix((new_data, (new_row, new_col)), shape=(3*rows,3*cols))
new_re = re.T @ re 
arr = re.toarray()
print(arr)
spchol(re.tocsc())
spchol(new_re.tocsc())