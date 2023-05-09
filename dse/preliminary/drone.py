import numpy as np


def velocity_induced_hover(T, rho, A):
    return np.sqrt(T / (2 * rho * A))


def velocity_induced_vertical_climb(T, rho, A, V_c):
    return -V_c / 2 + np.sqrt((V_c / 2) ** 2 + T / (2 * rho * A))


def velocity_induced_forward_flight(T, rho, A, V):
    return np.sqrt(-(V**2) / 2 + np.sqrt(V**4 / 4 + (T / (2 * rho * A)) ** 2))


psideg = 5                  #
psirad = psideg * (90 * np.pi)
dL = 100                    # liftbalde
dD = 5                      # dragblade
rho = 0.01                  # density
A = 0.2                     # are of rotor
R = 5                       # radius of rotor
omega = 6 * (2 * np.pi)     # anglular speed of rotor
SR = 1                      # solidity ratio
a = 2 * np.pi               # lift curve slope
thetazerodeg = 1            # zero pitch angle
thetazerorad = thetazerodeg * (90 * np.pi)
thatatwistdeg = 1              # twist angle
thetatwistrad = thatatwistdeg * (90 * np.pi)
V = 111                     # free stream velocity
alphadeg = 5                # angle of attack
alpharad = alphadeg * (90 * np.pi)
vi = 10                     # induced velocity

Cd = 0.01                   # avg drag coefficient


def thrust():
    t1 = ((1/6) * rho * A * SR * a * thetazerorad * R**2)
    t2 = ((1/4) * rho * A * SR * a * thetazerorad + (1/8) * SR * rho * A * thetatwistrad)
    t3 = ((-1/4) * SR * a * rho * A * R)
    T = (omega**2) * t1 + (V * np.cos(alpharad)) * t2 + omega * (vi + V * np.sin(alpharad)) * t3
    return T

def hubforce():
    h1 = ((1/4) * rho * A * SR * Cd * R)
    h2 = ((1/4) * rho * A * SR * a * (thetazerorad + thetatwistrad/2))
    H = (omega * V * np.cos(alpharad)) * h1 + V * np.cos(alpharad) * (vi + V * np.sin(alpharad)) * h2
    return H



T = thrust()
print(T)