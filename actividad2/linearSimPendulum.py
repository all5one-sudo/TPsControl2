import numpy as np
from control.matlab import *
import matplotlib.pyplot as plt

m = 0.1
F = 0.1
l = 1
g = 9.8
M = 1.5

# Matrices
A = np.array([[0, 1, 0, 0],
              [0, -F/M, -m*g/M, 0],
              [0, 0, 0, 1],
              [0, F/(l*M), g*(M+m)/(l*M), 0]])

B = np.array([[0],
              [1/M],
              [0],
              [-1/(l*M)]])

C = np.array([[1, 0, 0, 0],
              [0, 0, 1, 0]])

D = np.array([[0],
              [0]])

# Creación del modelo en espacio de estados
olSystem = ss(A, B, C, D)

# Polos a lazo abierto
olPoles = np.linalg.eigvals(A)

# Matriz de controlabilidad
M = ctrb(A, B)

# Variables para LQR
Q = np.dot(C.T, C)
R = 1
K, _, _ = lqr(A, B, Q, R)

# Nuevas matrices a implementar
Ac = A - np.dot(B, K)
Bc = B
Cc = C
Dc = D

# Sistema a lazo cerrado
clSys = ss(Ac, Bc, Cc, Dc)

# Tiempo de integración y simulación
h = 1e-4
simTime = 50
t = np.arange(0, simTime, h)

# Referencia del sistema
r = 10 * np.ones_like(t)

# Simulación
y, t, x = lsim(clSys, r, t)
delta = y[:, 0]
phi = y[:, 1]

# Gráfica
fig, ax = plt.subplots()
ax.plot(t, delta, label='Posición del carro [m]')
ax.set_xlabel('Tiempo [seg]')
ax.set_ylabel('Posición del carro [m]')
ax2 = ax.twinx()
ax2.plot(t, phi, 'r', label='Ángulo del péndulo [rad]')
ax2.set_ylabel('Ángulo del péndulo [rad]')
plt.title('Respuesta a entrada -10m con LQR')
plt.grid(True)
plt.show()
