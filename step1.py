import matplotlib.pyplot as plt
import numpy as np

xstart = 0
vstart = 0
tend = 50
dt = 0.01
x = xstart
v = vstart

m = 1
mu = 1
#1.38064852(79)×10−23
kb = 1
T = 300
D = mu * kb * T

gam = .19


L_v = []
L_x = []
L_t = []

def F(x):
    # https://en.wikipedia.org/wiki/Langevin_dynamics; R(t) is deze delta
    # wat is de standaarddev hier?
    # mu, sigma
    delta = np.random.normal(0, .5)
    fr = 2 * D * gam**2 * delta # delta ding

    F = x/2 - (x**3)/2 + fr
    return F

for t in np.arange(0, tend, dt):
    Fo = F(x)
    x = x + v * dt + Fo / (2 * m) * dt**2
    Fn = F(x)
    v = v + (Fo + Fn)/(2*m) * dt
    L_x.append(x)
    L_v.append(v)
    L_t.append(t)
    #print(x)

plt.plot(L_t,L_x)
#plt.plot(L_v, L_t)
plt.show()
