import matplotlib.pyplot as plt
import numpy as np

xstart = 0
vstart = 0
tend = 500
dt = 0.01
x = xstart
v = vstart
m = 27e-12

L_v = []
L_x = []
L_t = []

def F(x):
    return x/2 - (x**3)/2

for t in np.arange(0, tend, dt):
    Fo = F(x)
    x = x + v * dt + Fo / (2 * m) * dt**2
    Fn = F(x)
    v = v + (Fo + Fn)/(2*m) * dt
    L_x.append(x)
    L_v.append(v)
    L_t.append(t)

plt.plot(L_t,L_x)
#plt.plot(L_v, L_t)
plt.show()
