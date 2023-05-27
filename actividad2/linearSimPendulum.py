import numpy as np
from control.matlab import *
import matplotlib.pyplot as plt
from scipy.signal import place_poles, StateSpace, cont2discrete

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

# Gráfica del sistema real
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

# Planteo del observador
obsPoles = np.array([-6, -7, -8, -9])

# Matriz L
L = place_poles(A.T, C.T, obsPoles).gain_matrix.T

# Matrices de estado ampliadas
Aamp = np.block([[A - B @ K, B @ K], [np.zeros_like(A), A - L @ C]])
Bamp = np.block([[B], [np.zeros_like(B)]])
Camp = np.block([[Cc, np.zeros_like(Cc)]])
Damp = np.array([[0], [0]])

# Sistema a lazo cerrado
obsSysCL = StateSpace(Aamp, Bamp, Camp, Damp)

# Simulación lineal del observador
t1, y1, x1 = obsSysCL.output(U=r, T=t)

# Gráfica del sistema observado
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(t1, y1[:, 0])
ax2.plot(t1, y1[:, 1], 'r')
ax1.set_ylabel('Posición del carro observado [m]')
ax2.set_ylabel('Ángulo del péndulo observado [rad]')
ax1.set_xlabel('Tiempo [seg]')
plt.title('Respuesta a entrada -10m con LQR observada')
plt.grid(True)
plt.show()
