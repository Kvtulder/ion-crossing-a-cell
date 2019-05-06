# Define constants:

D = 1  # TODO: Einstein relation coefficient
GAMMA = 0.19
M = 27e-12  # weight of a cell, https://physics.aps.org/synopsis-for/10.1103/PhysRevLett.109.118105


# forces:

def F(x):
    return x / 2 - (x**3) / 2


def frictionalForce(x, t):
    # TODO implement calculation of force
    return 1;


def gradientForce(v):
    # TODO implement calculation of force
    return 1;


def environmentalInfluence():
    # TODO implement calculation of force
    return 1;

# calculates the average kinetic energy
def kineticEnergy():
    # TODO implement method
    return 1