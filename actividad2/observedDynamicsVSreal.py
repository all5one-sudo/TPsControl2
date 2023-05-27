import numpy as np
import matplotlib.pyplot as plt
from scipy.io import *

h = 1e-4
simTime = 50
t = np.arange(0,simTime,h)

deltaFile = loadmat('actividad2/deltaNoLineal.mat')
deltaPFile = loadmat('actividad2/deltaPNolineal.mat')
phiFile = loadmat('actividad2/phiNoLineal.mat')
phiPFile = loadmat('actividad2/phiPNoLineal.mat')
uNFile = loadmat('actividad2/uNoLineal.mat')

deltaOFile = loadmat('actividad2/deltaObs.mat')
deltaPOFile = loadmat('actividad2/deltaPObs.mat')
phiOFile = loadmat('actividad2/phiObs.mat')
phiOPFile = loadmat('actividad2/phiPO.mat')
uNOFile = loadmat('actividad2/uObs.mat')

delta = deltaFile['delta']
deltaP = deltaPFile['deltaP']
phi = phiFile['phi']
phiP = phiPFile['phiP']
u = uNFile['u']
deltaObs = deltaOFile['deltaO']
deltaPObs = deltaPOFile['deltaPO']
phiObs = phiOFile['phiO']
phiPObs = phiOPFile['phiPO']
uObs = uNOFile['u']


# Comparación del sistema real con el observado
'''plt.subplot(321)
plt.plot(t,np.squeeze(delta),label='Real')
plt.plot(t,np.squeeze(deltaObs),label='Observación')
plt.title('Desplazamiento del carro')
plt.ylabel('Distancia [m]')
plt.xlabel('Tiempo [seg]')
plt.legend(loc='best')
plt.grid()
plt.tight_layout()

plt.subplot(322)
plt.plot(t,np.squeeze(deltaP),label='Real')
plt.plot(t,np.squeeze(deltaPObs),label='Observación')
plt.title('Velocidad del carro')
plt.ylabel('Velocidad [m/s]')
plt.xlabel('Tiempo [seg]')
plt.legend(loc='best')
plt.grid()
plt.tight_layout()

plt.subplot(323)
plt.plot(t,np.squeeze(phi),label='Real')
plt.plot(t,np.squeeze(phiObs),label='Observación')
plt.title('Ángulo del carro')
plt.ylabel('Ángulo [rad]')
plt.xlabel('Tiempo [seg]')
plt.legend(loc='best')
plt.grid()
plt.tight_layout()

plt.subplot(324)
plt.plot(t,np.squeeze(phiP),label='Real')
plt.plot(t,np.squeeze(phiPObs),label='Observación')
plt.title('Velocidad angular')
plt.ylabel('Velocidad [rad/s]')
plt.xlabel('Tiempo [seg]')
plt.legend(loc='best')
plt.grid()
plt.tight_layout()

plt.subplot(313)
plt.plot(t,np.squeeze(u),label='Real')
plt.plot(t,np.squeeze(uObs),label='Observación')
plt.title('Acción de control')
plt.ylabel('u')
plt.xlabel('Tiempo [seg]')
plt.legend(loc='best')
plt.grid()
plt.tight_layout()'''

# Error de observación
plt.subplot(321)
plt.plot(t,(np.squeeze(delta)-np.squeeze(deltaObs)))
plt.title('Error en desplazamiento del carro')
plt.ylabel('Error')
plt.xlabel('Tiempo [seg]')
plt.grid()
plt.tight_layout()

plt.subplot(322)
plt.plot(t,(np.squeeze(deltaP)-np.squeeze(deltaPObs)))
plt.title('Error en velocidad del carro')
plt.ylabel('Error')
plt.xlabel('Tiempo [seg]')
plt.grid()
plt.tight_layout()

plt.subplot(323)
plt.plot(t,(np.squeeze(phi)-np.squeeze(phiObs)))
plt.title('Error en ángulo del carro')
plt.ylabel('Error')
plt.xlabel('Tiempo [seg]')
plt.grid()
plt.tight_layout()

plt.subplot(324)
plt.plot(t,(np.squeeze(phiP)-np.squeeze(phiPObs)))
plt.title('Error en velocidad angular')
plt.ylabel('Error')
plt.xlabel('Tiempo [seg]')
plt.grid()
plt.tight_layout()

plt.subplot(313)
plt.plot(t,(np.squeeze(u)-np.squeeze(uObs)))
plt.title('Error en acción de control')
plt.ylabel('Error')
plt.xlabel('Tiempo [seg]')
plt.grid()
plt.tight_layout()

plt.show()

