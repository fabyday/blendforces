import viewer 
import sys , os 


import numpy as np 



if __name__ == "__main__":
    app = viewer.QApplication(sys.argv)
    anim_data = np.load("./testanim_out/testanim.npy")
    data_path = "D:\\lab\\2022\\mycode\\FaceCaptureWithIK\\data\\ICT-data"
    neutral_pth = os.path.join(data_path, "generic_neutral_mesh.obj")
    win = viewer.MainWindow(neutral_pth, pipe=False, framerate=60)
    
    win.add_animation(anim_data)
    
    win.show()
    sys.exit(app.exec_())
    print("ex")