import numpy as np 

import mesh as mm 


import typing

class Blendshapes:
    def __init__(self, neutral : mm.Mesh, expression : typing.List[mm.Mesh]):
        self.__m_expression : list[mm.Mesh] = expression 
        self.__m_neutral  : mm.Mesh = neutral

    def build(self):
        self.__m_neutral.build_halfege()
        self.__B = np.hstack([ ( exp.v - self.__m_neutral.v).reshape(-1,1) for exp in self.__m_expression])
        
    def make_pose_by_weight(self, w):
        if not hasattr(self, "__B"):
            self.build()
        return self.__m_neutral.v + (self.__B @ w.reshape(-1,1)).reshape(-1, 3)

    def expression_pose(self, idx = None ):
        if idx is not None :
            return np.hstack([ (self.__m_neutral.v[idx, :] - exp.v[idx, :]).reshape(-1,1) for exp in self.__m_expression])
        return self.__B 
    
    def neutral_pose(self):
        return self.__m_neutral.v
    
    def neutral_mesh(self):
        return self.__m_neutral
        


