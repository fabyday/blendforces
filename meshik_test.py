import meshik as mk 
import subprocess


async def run_sub(queue, loop)    :
    print(loop, "ewew")
    proc = await loop.run_in_executor(
            None, 
            subprocess.Popen, 
            "python ./viewer.py", 0, None,
            subprocess.PIPE,       # stdin,
            
        )
    # proc = subprocess.Popen("python ./viewer.py", stdin=subprocess.PIPE)
    ii = 0 
    
    while True:
        data = await queue.get()
        proc.stdin.write(data)
        proc.stdin.flush()
        # print(f"put {ii}")
        ii+=1

def run_main(queue, loop):
    # Example usage
    import os, glob
    import mesh as mm
    meshik = mk.MeshIKConstraint()
    
    
    data_path = "D:\\lab\\2022\\mycode\\FaceCaptureWithIK\\data\\ICT-data"
    neutral_pth = os.path.join(data_path, "generic_neutral_mesh.obj")
    shapes_path = os.path.join(data_path, "shapes")
    file_pths = glob.glob(os.path.join(shapes_path, "**.obj"))
    neutral = mm.Mesh()
    neutral.load_from_file(neutral_pth)
    bs_list = []
    for fpth in file_pths:
        m = mm.Mesh()
        m.load_from_file(fpth)
        bs_list.append(m)
        
    lmk_idx = [1278,1272,12,1834,243,781,2199,1447,966,3661,4390,3022,2484,4036,2253,3490,3496,268,493,1914,2044,1401,3615,4240,4114,2734,2509,978,4527,4942,4857,1140,2075,1147,4269,3360,1507,1542,1537,1528,1518,1511,3742,3751,3756,3721,3725,3732,5708,5695,2081,0,4275,6200,6213,6346,6461,5518,5957,5841,5702,5711,5533,6216,6207,6470,5517,5966,]

    meshik.add_ref(neutral)
    meshik.add_examples(bs_list)
    meshik.precompute()
    import numpy as np 
    
    print("mesh ik loaded")
    print("mesh ik data load")
    datas = []
    for pth in glob.glob("./exported_objs/**.obj"):
        m = mm.Mesh()
        m.load_from_file(pth)
        datas.append(m.v[lmk_idx, :])
    print("mesh ik data load end")
    
    m_list = []
    R, C, RM  = meshik.getCMatrixAndReducedG(lmk_idx)
    ws = np.zeros(len(bs_list))
    for i in range(len(datas)):
        print("solving ", i)
        result = meshik.solve_nonlinear_from_lmk(RM, C, datas[i], lmk_idx)
        # result = meshik.TwLinear(ws)
        asyncio.run_coroutine_threadsafe(queue.put(result), loop)

        m_list.append(result)
    np.save("testanim_out/testik2.npy", np.array(m_list))
    m_list = np.load("testanim_out/testik2.npy")
    import viewer ,sys
    app = viewer.QApplication(sys.argv)

    win = viewer.MainWindow(neutral_pth, pipe=False, framerate=60)
    win.add_animation(m_list)
    
    win.show()
    sys.exit(app.exec_())
    print("ex")
    
if __name__ == "__main__":
    import threading
    import asyncio
    queue = asyncio.Queue()
    loop = asyncio.get_event_loop()
    print(loop, "init")
    calculation_thread = threading.Thread(target=run_main, args=(queue,loop))
    calculation_thread.daemon = True  # Ensure the thread exits when the main program exits
    calculation_thread.start()

    loop.run_until_complete(run_sub(queue,loop))


    import multiprocessing as mp 