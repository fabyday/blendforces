from PyQt5 import QtCore      # core Qt functionality
from PyQt5 import QtGui       # extends QtCore with GUI functionality
from PyQt5 import QtOpenGL    # provides QGLWidget, a special OpenGL QWidget

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import arcball
import gldeps
from OpenGL.WGL import *
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import OpenGL.GL as gl        # python wrapping of OpenGL
from OpenGL import GLU        # OpenGL Utility Library, extends OpenGL functionality
from OpenGL.arrays import vbo
import numpy as np
import sys                    # we'll need this later to run our Qt application
from PyQt5.QtCore import pyqtSignal, QObject
import time

class PipeThread(QThread):
    new_value = pyqtSignal(np.ndarray)
    def __init__(self , vsize, parent= None ):
        super().__init__(parent=parent)
        self.v_shape = vsize
    def run(self):
        read_byte = self.v_shape[0]*self.v_shape[1]*np.dtype(np.float64).itemsize
        while True :
            x = sys.stdin.buffer.read(read_byte)
            x = np.frombuffer(x, dtype=np.float64).reshape(-1, 3)
            self.new_value.emit(x)
            self.msleep(500)
import igl 

class GLWidget(QtOpenGL.QGLWidget):
    
    
    
    def __init__(self, neutral_mesh_path  ,parent=None):
        self.parent = parent
        QtOpenGL.QGLWidget.__init__(self, parent)
        self.v, self.f  = igl.read_triangle_mesh(neutral_mesh_path)


        self.setMouseTracking(True)
        self.zoom = 0
        self.setRotX(0)
        self.setRotY(0)
        self.setRotZ(0)


    @pyqtSlot(np.ndarray)
    def update(self, x):
        self.vert_vbo.set_array(x) 
        self.v = x
        self.update_normal()


    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(0, 0, 0))    # initialize the screen to blue
        gl.glEnable(gl.GL_DEPTH_TEST)                  # enable depth testing
        glFrontFace(GL_CCW)
        glEnable(GL_CULL_FACE)
        glEnable(GL_LIGHTING)
        self.ambientLight = np.array([0.5,0.5,0.5,1.0])
        self.diffuseLight = np.array([0.5, 0.5 ,0.5 ,1.0])
        self.specular = np.array([1.0, 1.0, 1.0, 1.0])
        self.lightPosition = np.array([0.0,0.0, 1.0, 0.0])
        glLightfv(GL_LIGHT0, GL_AMBIENT, self.ambientLight)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, self.diffuseLight)
        glLightfv(GL_LIGHT0, GL_SPECULAR, self.specular)
        glLightfv(GL_LIGHT0, GL_POSITION, self.lightPosition)
        glEnable(GL_LIGHT0)

        self.initGeometry()
  


    def setRotX(self, val):
        self.rotX = val

    def setRotY(self, val):
        self.rotY = val

    def setRotZ(self, val):
        self.rotZ = val



    def update_normal(self):
        v1 = self.v[self.f[:, 0], :]
        v2 = self.v[self.f[:, 1], :]
        v3 = self.v[self.f[:, 2], :]

        v12 = v2 -v1
        v13 = v3 -v1 
        f_normals = np.cross(v12, v13)
        v_n = np.zeros_like(self.v)
        v_n_denom = np.zeros((len(self.v),1))
        for fi, fn in enumerate(f_normals):
            v1, v2, v3 = self.f[fi]
            v_n[v1] += fn
            v_n[v2] += fn
            v_n[v3] += fn
            v_n_denom[v1] += 1
            v_n_denom[v2] += 1
            v_n_denom[v3] += 1
        self.v_n = v_n/v_n_denom
        self.v_n /= np.linalg.norm(self.v_n, axis=-1).reshape(-1,1)




        

    
    def paintGL(self):
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
        mean = np.mean(self.v, axis=0)
        bmax = np.max(self.v, axis=0)
        bmin = np.min(self.v, axis=0)
        self.center = (bmax + bmin)*0.5
        self.scale = np.linalg.norm(bmax - self.center)
        gl.glPushMatrix()  
        glLoadIdentity();   # push the current matrix to the current stack
        x,y,z = (self.center).ravel()
        
        
        gluLookAt(x,y,z+ 100 - self.zoom, x,y,z,0.0 ,1.0 ,0)
        gl.glTranslate(0.0, 0.0, -10.0)    # third, translate cube to specified depth
        # gl.glScale(self.scale, self.scale, self.scale)
        gl.glRotatef(self.rotX,1,0,0)
        gl.glRotatef(self.rotY, 0, 1, 0)
        gl.glRotatef(self.rotZ, 0,0,1)
        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
        gl.glEnableClientState(gl.GL_COLOR_ARRAY)

        gl.glVertexPointer(3, gl.GL_FLOAT, 0, self.v.reshape(1,-1).astype(np.float32))
        c = np.ones_like(self.v)
        gl.glNormalPointer(gl.GL_FLOAT, 0, self.v_n.reshape(1,-1).astype(np.float32))
        gl.glColorPointer(3, gl.GL_FLOAT, 0, c.reshape(1,-1).astype(np.float32))
        gl.glDrawElements(gl.GL_TRIANGLES, len(self.f)*3, gl.GL_UNSIGNED_INT, self.f.reshape(1,-1).astype(np.uint))

        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
        gl.glDisableClientState(gl.GL_NORMAL_ARRAY)
        

        gl.glPopMatrix()    # restore the previous modelview matrix

    def resizeGL(self, width, height):
        gl.glViewport(0, 0, width, height)
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        aspect = width / float(height)

        GLU.gluPerspective(45.0, aspect, 1.0, 1000.0)
        gl.glMatrixMode(gl.GL_MODELVIEW)
    
    def set_data(self, v,f):
        self.v = v 
        self.f = f


    def wheelEvent(self,event):
        self.zoom += event.angleDelta().y()/30


    def initGeometry(self):
        
        self.shader = gldeps.Shader()
        self.shader.compile(gldeps.v_shader_src, gldeps.p_shader_src)
        
        self.vert_vbo = vbo.VBO(self.v.reshape(1,-1).astype(np.float32))
        self.f_vbo = vbo.VBO(self.f.astype(np.uint))
        self.update_normal()

        


class MainWindow(QMainWindow):

    def __init__(self, neutral_mesh_path):
        QMainWindow.__init__(self)    # call the init for the parent class
        v, f = igl.read_triangle_mesh(neutral_mesh_path)
        shape = v.shape
        self.resize(800, 800)
        self.setWindowTitle('viewer')
        self.worker_thread = PipeThread(shape, self)
        self.worker_thread.start()

        self.glWidget = GLWidget(neutral_mesh_path, self)
        self.initGUI()
        
        timer = QtCore.QTimer(self)
        timer.setInterval(20)   # period, in milliseconds
        timer.timeout.connect(self.glWidget.updateGL)
        timer.start()
        self.worker_thread.new_value.connect(self.glWidget.update)

    def initGUI(self):
        central_widget = QWidget()
        gui_layout = QVBoxLayout()
        central_widget.setLayout(gui_layout)

        self.setCentralWidget(central_widget)

        gui_layout.addWidget(self.glWidget)

        sliderX = QSlider(QtCore.Qt.Horizontal)
        sliderX.valueChanged.connect(lambda val: self.glWidget.setRotX(val))

        sliderY = QSlider(QtCore.Qt.Horizontal)
        sliderY.valueChanged.connect(lambda val: self.glWidget.setRotY(val))

        sliderZ = QSlider(QtCore.Qt.Horizontal)
        sliderZ.valueChanged.connect(lambda val: self.glWidget.setRotZ(val))
        
        gui_layout.addWidget(sliderX)
        gui_layout.addWidget(sliderY)
        gui_layout.addWidget(sliderZ)

if __name__ == '__main__':

    app = QApplication(sys.argv)
    
    data_path = "D:\\lab\\2022\\mycode\\FaceCaptureWithIK\\data\\ICT-data"
    neutral_pth = os.path.join(data_path, "generic_neutral_mesh.obj")
    win = MainWindow(neutral_pth)
    win.show()
    sys.exit(app.exec_())
    print("ex")