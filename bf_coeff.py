import numpy as np 
import igl 
import os 
import os.path as osp 


class CoeffFactory():
    def __init__(self, default_pth = None ):
        self.init = False 
        self.mesh_pth = default_pth if default_pth is not None else "neutral_mesh.obj"
    
    
    def _compute(self, save_pth):
        pass 
    
    def load_coeffs(self, chache_dir):
        """
        Load coefficients from a cache file.
        """
        if not osp.exists(chache_dir):
            print(f"Cache directory {chache_dir} does not exist. Computing coefficients...")
            self._compute(chache_dir)
        
        coeffs = np.load(chache_dir, allow_pickle=True).item()
        
        self.coeffs = coeffs
        self.init = True
        return coeffs
    
    
    
if __name__ == "__main__": 
    CoeffFactory().load_coeffs("./coeffs")
    