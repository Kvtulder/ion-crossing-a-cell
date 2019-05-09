import matplotlib.pyplot as plt
import numpy as np
import statistics
import math

xstart = 0
vstart = 0
tend = 1500

# uit t paper
dt = 0.01


m = 1
mu = 1
#m = 1.38064852 * 10**(-23)
kb = 1
T = 310
D = mu * kb * T

gam = .19

# constants for random force:
omega = 0.833
f0 = 0.16



def F(x,t,v):
    # https://en.wikipedia.org/wiki/Langevin_dynamics; R(t) is deze delta
    # wat is de standaarddev hier?
    # mu, sigma
    delta = np.random.normal(0,1)
    fr = 2 * D * gam**2 * delta # delta ding

    # dit is uit het paper
    #fr = f0 * math.cos(omega * t)

    # kasper; dit is de wrijving die we waren vergeten (nog niet hadden)
    friction = - gam * v

    F = x/2 - (x**3)/2 + fr  + friction
    return F

def V(x):
    return - (x ** 2) / 4 + (x ** 4) / 8

def simulate(plot = True):

    x = xstart
    v = vstart

    L_v = []
    L_x = []
    L_t = []

    Ekintot = 0

    for t in np.arange(0, tend, dt):
        Fo = F(x,t,v)
        x = x + v * dt + Fo / (2 * m) * dt**2
        Fn = F(x,t,v)
        v = v + (Fo + Fn)/(2*m) * dt

        Ekin = .5 * m * v**2

        Ekintot = Ekintot + Ekin

        if t > 1000:
            L_x.append(x)
            L_v.append(v)
            L_t.append(t)
            #print("hi")

    open = False
    startopen = 0
    periods = []

    #print(L_x)

    # calculates how long the channel stays open
    for i in range(0, len(L_x) - 1):
        if L_x[i] > 0 and not open:
            # channel just opened!
            open = True
            startopen = i
        elif L_x[i] < 0 and open:
            # channel just closed
            open = False
            # amount of frames between the two events times dt is out period
            periods.append((i - startopen) * dt)
            #print(i)
    # only plot if wanted
    if plot:
        # normal plot
        plt.subplot(2,1,1)
        plt.plot(L_t, L_x)
        plt.xlabel('tijd')
        plt.ylabel('plaats')
        plt.text(1,1, "average open time %d" % statistics.mean(periods))
        plt.subplot(2,1,1)
        plt.subplot(2,1,2)
        plt.plot(L_x, L_v)
        plt.plot(x = 0)
        plt.xlabel('plaats')
        plt.ylabel('snelheid')
        plt.show()

    print('Average kin energy',Ekintot/(tend/dt))

    return statistics.mean(periods)

periods = []

for time in range(3):
    print("Calculation + %d / 100" % time)
    periods.append(simulate(False))

print("average open time %f" %statistics.mean(periods))
simulate(True)
