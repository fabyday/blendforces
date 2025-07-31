import mesh 

import geometry_helper as geo 

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
        
        
    def add_examples(self, examples : list[mesh.Mesh]):
        self.m_examples = examples

    
    def precompute(self):
        if not hasattr(self, "m_ref"):
            raise  ValueError("Reference mesh is not set.")
        self.vn = geo.compute_vertex_normals(self.m_ref.v, self.m_ref.f)
        self.extended_f = np.concatenate([self.m_ref.f, len(self.m_ref.v)+np.arange(len(self.vn)) ],-1)
        
        
        
        
        
    def GMatrix(self):
        return 
    
 
    
    
    