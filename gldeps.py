import sys
import array
import ctypes
import struct
from typing import Any
from OpenGL.GL import *
from OpenGL.GLUT import *
v_shader_src = """
            #version 330 core
            layout (location = 0) in vec3 aPos;
            layout (location = 1) in vec3 anormal;

            uniform mat4 Q;
            uniform mat4 Rt;

            out vec3 Normal;
            out vec3 FragPos;

            void main()
            {
                gl_Position = Q*Rt*vec4(aPos, 1.0);
                //FragPos = vec3(Rt*vec4(aPos,1.0));
                //Normal = mat3(transpose(inverse(Rt)))*anormal;
                FragPos =aPos;
                Normal = anormal;
            }
            """

p_shader_src ="""
    #version 330 core
    out vec4 FragColor;

    uniform vec3 lightColor;
    uniform vec3 objectColor;
    uniform vec3 lightPos;

    in vec3 Normal;
    in vec3 FragPos;


    void main()
    {
        float ambientStrength = 0.1;
        vec3 ambient = ambientStrength * (lightColor);
        vec3 norm = normalize(Normal);
        vec3 lightDir = normalize(lightPos - FragPos);
        
        float diff = max(dot(norm, lightDir), 0.0);
        vec3 diffuse = diff*lightColor;
        vec3 result = (ambient + diffuse)*objectColor;
        FragColor = vec4(result, 0.0);
    }
    """
def load_shader(src: str, shader_type: int) -> int:
    shader = glCreateShader(shader_type)
    glShaderSource(shader, src)
    glCompileShader(shader)
    error = glGetShaderiv(shader, GL_COMPILE_STATUS)
    if error != GL_TRUE:
        info = glGetShaderInfoLog(shader)
        glDeleteShader(shader)
        raise Exception(info)
    return shader


class Shader:
    def __init__(self) -> None:
        self.program = glCreateProgram()

    def __del__(self) -> None:
        glDeleteProgram(self.program)

    def compile(self, vs_src: str =v_shader_src, fs_src: str = p_shader_src) -> None:
        vs = load_shader(vs_src, GL_VERTEX_SHADER)
        if not vs:
            return
        fs = load_shader(fs_src, GL_FRAGMENT_SHADER)
        if not fs:
            return
        glAttachShader(self.program, vs)
        glAttachShader(self.program, fs)
        glLinkProgram(self.program)
        error = glGetProgramiv(self.program, GL_LINK_STATUS)
        glDeleteShader(vs)
        glDeleteShader(fs)
        if error != GL_TRUE:
            info = glGetShaderInfoLog(self.program)
            raise Exception(info)

    def use(self):
        glUseProgram(self.program)

    def unuse(self):
        glUseProgram(0)




