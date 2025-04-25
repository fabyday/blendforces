import numpy as np
import matplotlib.pyplot as plt

# Mass class
class Mass:
    def __init__(self, pos, mass=1.0, fixed=False):
        self.pos = np.array(pos, dtype=np.float64)
        self.prev_pos = np.array(pos, dtype=np.float64)
        self.mass = mass
        self.fixed = fixed
        self.force = np.zeros(2)

    def apply_force(self, f):
        if not self.fixed:
            self.force += f

    def integrate(self, dt):
        if self.fixed:
            return
        acc = self.force / self.mass
        next_pos = self.pos + (self.pos - self.prev_pos) + acc * dt**2
        self.prev_pos = self.pos.copy()
        self.pos = next_pos
        self.force[:] = 0

# Spring class
class Spring:
    def __init__(self, m1, m2, rest_length, stiffness=1.0):
        self.m1 = m1
        self.m2 = m2
        self.rest_length = rest_length
        self.k = stiffness

    def apply_force(self):
        delta = self.m2.pos - self.m1.pos
        dist = np.linalg.norm(delta)
        if dist == 0:
            return
        force = self.k * (dist - self.rest_length) * (delta / dist)
        self.m1.apply_force(force)
        self.m2.apply_force(-force)

# 시뮬레이션 설정
m1 = Mass([0.0, 0.0], fixed=True)
m2 = Mass([1.0, 0.0])
spring = Spring(m1, m2, rest_length=1.0, stiffness=10.0)

masses = [m1, m2]
springs = [spring]

# 시뮬레이션 루프
dt = 0.01
positions = []

for _ in range(500):
    # Gravity
    for m in masses:
        m.apply_force(np.array([0, -9.8 * m.mass]))

    # Spring forces
    for s in springs:
        s.apply_force()

    # Integrate
    for m in masses:
        m.integrate(dt)

    positions.append(m2.pos.copy())

# 결과 시각화
positions = np.array(positions)
plt.plot(positions[:, 0], positions[:, 1])
plt.title('Mass-Spring Simulation')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.show()