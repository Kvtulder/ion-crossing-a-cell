import matplotlib.pyplot as plt
import numpy as np
import statistics
import math

xstart = 0
vstart = 0
tend = 1500

# uit t paper
dt = 0.01
#m = 27**(-20)
m = 1
#m = 1.38064852 * 10**(-23)
#kb = 1.38065*10**(-23)
kb = 1
T = 310
#T = 600
E_equ = .5 * kb * T
print(E_equ, ":equip")


gam = .19

# constants for random force:
omega = 0.833
f0 = 0.16

def F(x,t,v):
    # https://en.wikipedia.org/wiki/Langevin_dynamics; R(t) is deze delta
    # wat is de standaarddev hier?
    # mu, sigma
    delta = np.random.normal(0,1)
    if math.fabs(V(x))> 0.05:
        #print(x)
        #print(V(x), "pot")
        mu = v/V(x)
        mu = 1
        #print(mu)
        #mu = 1
    else:
        mu = 1
    D = mu * kb * T
    #fr = 2 * D * gam**2 * delta

    # dit is uit het paper
    fr = f0 * math.cos(omega * t)

    friction = - gam * v

    F = x/2 - (x**3)/2 + fr  + friction
    #print(F, "powerrrs")
    return F

def V(x):
    return - (x ** 2) / 4 + (x ** 4) / 8

def simulate(plot = True):

    x = xstart
    v = vstart

    L_v = []
    L_fabsv = []
    L_x = []
    L_t = []

    L_E = []
    L_Ekin = []
    L_Epot = []

    L_maxie= []
    L_vmaxie = []

    Ekintot = 0


    for t in np.arange(0, tend, dt):
        Fo = F(x,t,v)
        x = x + v * dt + Fo / (2 * m) * dt**2
        Fn = F(x,t,v)
        v = v + (Fo + Fn)/(2*m) * dt

        if t > 1000:
            L_x.append(x)
            L_v.append(v)
            L_fabsv.append(math.fabs(v))
            L_t.append(t)

            Ekin = .5 * m * v**2
            Epot = V(x)
            L_E.append(Ekin+Epot)
            L_Ekin.append(Ekin)
            L_Epot.append(Epot)

            Ekintot = Ekintot + Ekin

    open = False
    startopen = 0
    periods = []

    #maxwell boltzman
    vmin = float(min(L_v))
    vmax = float(max(L_v))
    Cmax = (m*10**2/(2*math.pi*kb*T))**.5
    for v in np.arange(vmin,vmax, 0.01):
        maxie = Cmax* math.exp(-m*10**2*v**2/(2*kb*T))
        L_maxie.append(maxie)
        L_vmaxie.append(v)

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
        plt.xlabel('time [s]')
        plt.ylabel('place [a.u.]')
        plt.text(1,1, "average open time %d" % statistics.mean(periods))
        plt.title("place-time plot")
        #plt.subplot(2,1,1)
        plt.subplot(2,1,2)
        plt.plot(L_x, L_v)
        plt.plot(x = 0)
        plt.xlabel('place [a.u.]')
        plt.ylabel('speed [a.u.]')

        plt.subplots_adjust(hspace=0.5)
        plt.title("phase-space plot")
        plt.show()

    #potentiaal plotje
    xmin = min(L_x)
    xmax = max(L_x)
    L_xpotentiaal= []
    L_Epot = []
    for x in np.arange(xmin,-xmin, 0.01):
        L_xpotentiaal.append(x)
        L_Epot.append(V(x))


    print('Average kin energy',Ekintot/((tend-1000)/dt))
    #plt.plot(L_t, L_E)
    #plt.plot(L_t, L_Ekin)
    #plt.plot(L_xpotentiaal, L_Epot)
    #plt.xlabel("x [a.u.]")
    #plt.ylabel("Epot [a.u.]")
    #plt.title("Potential of an ion motion in or out a cell through a protein channel")

    #plt.xlim(1000,1100)
    plt.show()
    #print(L_E)

    #plt.hist(L_v, bins = 100, normed = 1)
    #plt.plot(L_vmaxie, L_maxie)
    #plt.show()

    return statistics.mean(periods)

periods = []

for time in range(3):
    print("Calculation + %d / 100" % time)
    periods.append(simulate(False))

print("average open time %f" %statistics.mean(periods))
simulate(True)
