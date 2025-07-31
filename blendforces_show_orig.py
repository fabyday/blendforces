import viewer 
import sys , os , glob
import mesh as mm


import numpy as np 



if __name__ == "__main__":
    lmk_idx = [1278,1272,12,1834,243,781,2199,1447,966,3661,4390,3022,2484,4036,2253,3490,3496,268,493,1914,2044,1401,3615,4240,4114,2734,2509,978,4527,4942,4857,1140,2075,1147,4269,3360,1507,1542,1537,1528,1518,1511,3742,3751,3756,3721,3725,3732,5708,5695,2081,0,4275,6200,6213,6346,6461,5518,5957,5841,5702,5711,5533,6216,6207,6470,5517,5966,]
    datas = []
    for pth in glob.glob("./exported_objs/**.obj"):
        m = mm.Mesh()
        m.load_from_file(pth)
        datas.append(m.v)
    datas = np.array(datas)
    app = viewer.QApplication(sys.argv)
    data_path = "D:\\lab\\2022\\mycode\\FaceCaptureWithIK\\data\\ICT-data"
    neutral_pth = os.path.join(data_path, "generic_neutral_mesh.obj")
    
    win = viewer.MainWindow(neutral_pth, pipe=False, framerate=60)
    
    win.add_animation(datas)
    
    win.show()
    sys.exit(app.exec_())
    print("ex")