import scipy.sparse as sp 
import numpy as np 


data = np.arange(1, 10).reshape(3,3)
s = sp.coo_matrix(data, shape=(3,3))
d = s.data 
row = s.row
col = s.col
ar = np.arange(1,10).reshape(3,3)
print(data)
print(data@ar)
new_row = [] 
new_col = [] 
new_data = []
for r, c, dd in zip(row, col, d):
    for axis in range(3):
        new_row.append(3 * r + axis)
        new_col.append(3 * c + axis)
        new_data.append(dd)
re = sp.coo_matrix((new_data, (new_row, new_col)), shape=(9,9))
arr = re.toarray()
print((arr@ar.reshape(-1,1)).reshape(3,3))
print(arr)

    