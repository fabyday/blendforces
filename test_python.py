import numpy as np 



A1 = np.array([[1],[-1]])
A2 = np.array([[-1],[1]])
A1tA1 = A1 @ A1.T
Af1 = np.array([[1, 0, 0, -1, 0, 0],
               [0,1,0, 0,-1,0],
               [0,0,1, 0,0,-1],
               
               ])
Af2 = np.array([[-1,0,0, 1,0,0],
               [0,-1,0, 0,1,0],
               [0,0,-1, 0,0,1],])
print(Af1.T@Af1 + Af2.T@Af2)

Sf = np.array([[1, 0, 0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0]])
Sf2 = np.array([[0,0,0,1, 0, 0],[0,0,0,0,1,0,],[0,0,0,0,0,1]])
print("SS\n", Af1.T@Sf + Af2.T@Sf2)

print(np.kron(A1tA1, np.identity(3)))
A2tA2 = A2@A2.T 
s = np.kron(
(A1tA1 + A2tA2), np.identity(3))
print(s)


S1 = np.array([[1],[0]])

S2 = np.array([[0],[1]])

b= np.kron((A1 @ S1.T + A2@S2.T), np.identity(3))
print(b)


np.kron((A1 @ S1.T + A2@S2.T), np.identity(3))




