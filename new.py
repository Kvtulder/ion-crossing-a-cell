import numpy as np
import matplotlib.pyplot as plt


# define constants
x_start, v_start = 0, 0
gamma = 0.19
amu = 1.660539 * 10**-27
m = 39.0983 * amu
m = 1
T = 310
kb = 1.38064852 * 10**-23
kb = 3*10**3
D = kb * T / (gamma * m)
delta_t = 0.01

# we use a normal distribution with unit variance, so we use mean 0 and variance 1
mean = 0
variance = 1


# calculate the potential
def V(x_var):
    return -(x_var**2 / 4.0) + (x_var**4 / 8)


def x(x_prev, v_prev, random_1, random_2):
    return x_prev + delta_t * v_prev + A(x_prev, v_prev, random_1, random_2)


def v(x, x_prev, v_prev, random_1, random_2):
    return v_prev + (0.5 * delta_t) * (F_p(x) + F_p(x_prev)) \
           - (gamma / m) * A(x_prev, v_prev, random_1, random_2) \
           + (2 * D * gamma ** 2 * delta_t) * random_1

# F_p is calculated by deriving V to x * -1
def F_p(x_var):
    return x_var/2 - (x_var**3)/2


def A(x, v, random_1, random_2):
    p1 = (0.5 * delta_t ** 2) * (F_p(x) - (gamma / m) * v)
    p2 = ((2 * D * (gamma ** 2)) ** 0.5) * delta_t ** 1.5
    p3 = 0.5 * random_1 + 0.5 * ((3)**(-0.5)) * random_2

    return p1 + p2 * p3

# returns two random variables using the standard normal distribution
def normal(mean, variance):
    s = np.random.normal(mean, variance, 2)
    return s[0], s[1]

t_list = [0]
x_list = [x_start]
v_list = [v_start]

for t in np.arange(0, 50, delta_t):

    random_1, random_2 = normal(mean, variance)

    x_prev = x_list[-1]
    v_prev = v_list[-1]

    x_new = x(x_prev, v_prev, random_1, random_2)
    v_new = v(x_new, x_prev, v_prev, random_1, random_2)

    x_list.append(x_new)
    v_list.append(v_new)
    t_list.append(t)

plt.plot(t_list, x_list)
plt.show()
