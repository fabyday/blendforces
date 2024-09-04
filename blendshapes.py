import numpy as np 

import mesh as mm 




class Blendshapes:
    def __init__(self, neutral : mm.Mesh, expression : list[mm.Mesh]):
        self.__m_expression : list[mm.Mesh] = expression 
        self.__m_neutral  : mm.Mesh = neutral


    def build(self):
        self.__B = np.hstack([ self.__m_neutral.v - exp.v for exp in self.__m_expression])
        
    def make_pose_by_weight(self, w):
        return self.__m_neutral.v + self._B @ w.reshape(-1,1)

    def expression_pose(self):
        return self.__B 
    
    def neutral_pose(self):
        return self.__m_neutral.v
        


