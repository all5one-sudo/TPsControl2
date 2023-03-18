from control.matlab import *
import numpy as np
import matplotlib.pyplot as plt

num = [1]
den = [1,1]
sys = tf(num,den)

t = np.arange(0,20,0.01)
w0 = 3
u = np.cos(w0*t)
y,t,x = lsim(sys, u, t)

plt.plot(t,y)
plt.grid()
plt.show()