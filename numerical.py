# this file is used for all the integration techniques
import physics


def verletX(x_prev, v_prev, dt):

    return x_prev + v_prev * dt + (physics.F(x) / physics.M) * dt ** 2

def verletV(v_prev):
    return v_prev + (physics.F())