import matplotlib.pyplot as plt
import numpy as np
import statistics

xstart = 0
vstart = 0
tend = 500
dt = 0.01


m = 1
mu = 1
#1.38064852(79)×10−23
kb = 1
T = 300
D = mu * kb * T

gam = .19




def F(x):
    # https://en.wikipedia.org/wiki/Langevin_dynamics; R(t) is deze delta
    # wat is de standaarddev hier?
    # mu, sigma
    delta = np.random.normal(0, .5)
    fr = 2 * D * gam**2 * delta # delta ding

    F = x/2 - (x**3)/2 + fr
    return F

def V(x):
    return - (x ** 2) / 4 + (x ** 4) / 8

def simulate(plot = True):

    x = xstart
    v = vstart

    L_v = []
    L_x = []
    L_t = []

    for t in np.arange(0, tend, dt):
        Fo = F(x)
        x = x + v * dt + Fo / (2 * m) * dt**2
        Fn = F(x)
        v = v + (Fo + Fn)/(2*m) * dt
        L_x.append(x)
        L_v.append(v)
        L_t.append(t)

    open = False
    startopen = 0
    periods = []

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

    # only plot if wanted
    if plot:
        # normal plot
        plt.subplot(2,1,1)
        plt.plot(L_t, L_x)
        plt.text(1,1, "average open time %d" % statistics.mean(periods))
        plt.subplot(2,1,1)
        plt.subplot(2,1,2)
        plt.plot(L_x, L_v)
        plt.show()

    return statistics.mean(periods)


periods = []

for time in range(3):
    print("Calculation + %d / 100" % time)
    periods.append(simulate(False))

print("average open time", statistics.mean(periods))
simulate(True)

