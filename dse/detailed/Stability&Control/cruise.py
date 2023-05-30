import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# atmospheric conditions
rho = 0.02
m = 3000
g = 3.71
W0 = m * g
I = [81250, 19750, 91500]  # mass moments of inertia ixx, iyy, izz
Vcruise = 400/3.6


# plane design parameters
downwash = 3        # in deg
instalation = 0     # in deg
area = [56, 56, 15] # left wing, right wing, stabilizer

def definegeometry():
    '''

    Returns: the position of each force with respect to the aircraft body reference frame

    '''
    Rwingleft = np.array([1, -10, -1])
    Rwingright = np.array([1, 10, -1])
    Rstabilizer = np.array([-10, 0, -2])
    Raileronleft = np.array([0, -16, -1])
    Raileronright = np.array([0, 16, -1])
    Rrudder = np.array([-10, 0, -3])
    Cg = np.array([0, 0, 0])
    Rthrustleft = np.array([2, -17, -1])
    Rthrustright = np.array([2, 17, -1])
    R = np.array([Rwingleft, Rwingright, Rstabilizer, Raileronleft, Raileronright, Rrudder, Rthrustright, Rthrustleft, Cg])
    return R

def getcoefficients(aoal, aoar, aoah):
    '''

    Args:
        aoal: angle of attack left wing
        aoar: angle of attack right wing
        aoah: angle of attack horizontal stabilizer

    Returns:
        clwl: lift coefficient of left wing
        clwr: lift coefficient of right wing
        clh: lift coefficient of horizontal stabilizer
        cdwl: drag coefficient of left wing
        cdwr: drag coefficient of right wing
        cdh: drag coefficient of horizontal stabilizer
    '''
    aoalistw = np.arange(-5, 15, 0.1)
    aoalisth = np.arange(-5, 15, 0.1)
    cllistw = np.radians(aoalistw) * 2 * np.pi   # will later be a csv file
    cllisth = np.radians(aoalisth) * 2 * np.pi   # will later be a csv file
    cldivcdw = np.ones(len(aoalistw)) * 20   # will later be a csv file
    cldivcdh = np.ones(len(aoalisth)) * 20   # will later be a csv file
    clwl = np.interp(aoal, aoalistw, cllistw)
    clwr = np.interp(aoar, aoalistw, cllistw)
    clh = np.interp(aoah, aoalisth, cllisth)
    cdwl = clwl / np.interp(aoal, aoalistw, cldivcdw)
    cdwr = clwr / np.interp(aoar, aoalistw, cldivcdw)
    cdh = clh / np.interp(aoah, aoalisth, cldivcdh)
    return clwl, clwr, clh, cdwl, cdwr, cdh


def aerodynamic(aoa, v, beta):
    '''

    Args:
        aoa: angle of attack of the left wing, right wing and horizontal stabilizer
        v: speed of left wing, right wing and horizontal stabilizer
        beta: sideslip angle

    Returns:
        wl: force vector left wing
        wr: force vector right wing
        h: force vector horizontal stabilizer
    '''
    # aoa includes the effective angle of attack due to the pitch and angular velocity
    clwl, clwr, clh, cdwl, cdwr, cdh = getcoefficients(aoa[0], aoa[1], aoa[2])
    Lwl = 0.5 * rho * area[0] * clwl * v[0]**2
    Lwr = 0.5 * rho * area[1] * clwr * v[1]**2
    Lh = 0.5 * rho * area[2] * clh * v[2]**2
    Dwl = 0.5 * rho * area[0] * cdwl * v[0]**2
    Dwr = 0.5 * rho * area[1] * cdwr * v[1]**2
    Dh = 0.5 * rho * area[2] * cdh * v[2]**2

    wl = np.array([-Dwl, 0, -Lwl])
    wr = np.array([-Dwr, 0, -Lwr])
    h = np.array([-Dh, 0, -Lh])
    wl = np.matmul(aero_to_body(aoa[3], beta), wl).T
    wr = np.matmul(aero_to_body(aoa[3], beta), wr).T
    h = np.matmul(aero_to_body(aoa[3], beta), h).T
    return wl, wr, h

def thrust(v, aoa, a, Tl, Tr):
    if a[0]>0.05:
        #accelerating too fast
        Tl = Tl - 50
        Tr = Tr - 50
    elif a[0]<-0.05:
        #deaccelerating to fast
        Tr = Tr + 50
        Tl = Tl + 50
    elif v[0] > Vcruise * np.cos(np.radians(aoa[3])) and a[0] > 0:
        #accelerating away from cruise
        Tr = Tr - 25
        Tl = Tl - 25
    elif v[0] < Vcruise * np.cos(np.radians(aoa[3])) and a[0] < 0:
        #deaccelerating away from cruise
        Tr = Tr + 25
        Tl = Tl + 25
    elif v[0] < Vcruise * np.cos(np.radians(aoa[3])) and a[0] > 0:
        #accelerating towards cruise
        Tr = Tr - 25 * abs(v[0] - Vcruise * np.cos(np.radians(aoa[3])))/abs(Vcruise * np.cos(np.radians(aoa[3])))
        Tl = Tl - 25 * abs(v[0] - Vcruise * np.cos(np.radians(aoa[3])))/abs(Vcruise * np.cos(np.radians(aoa[3])))
    elif v[0] > Vcruise * np.cos(np.radians(aoa[3])) and a[0] < 0:
        #deacelerating towards cruise
        Tr = Tr + 25 * abs(v[0] - Vcruise * np.cos(np.radians(aoa[3]))) / abs(Vcruise * np.cos(np.radians(aoa[3])))
        Tl = Tl + 25 * abs(v[0] - Vcruise * np.cos(np.radians(aoa[3]))) / abs(Vcruise * np.cos(np.radians(aoa[3])))
    else:
        Tl = Tl
        Tr = Tr

    return Tl, Tr

def control():
    return 1

def moments(R, F):
    '''

    Args:
        R: list of the distances for left wing, right wing, horizontal stabilizer, aileron left, aileron right, rudder and cg
        F: list of forces for the same

    Returns:
        M: sum of moment arount the X axis, Y axis and Z axis
    '''

    M = np.zeros((len(R), 3))
    for i in range(len(R)):
        M[i] = np.cross(R[i], F[i])
    M = np.array([sum(M.T[0]), sum(M.T[1]), sum(M.T[2])])
    return M

def accelerations(F, M, I):

    a = [sum(F.T[0])/m, sum(F.T[1])/m, sum(F.T[2])/m]
    ang = [M[0]/I[0], M[1]/I[1], M[2]/I[2]]
    return a, ang



def aero_to_body(a, b):
    a = np.radians(a)
    b = np.radians(b)
    T = np.array([[np.cos(b) * np.cos(a), np.sin(b), np.cos(b) * np.sin(a)],
                  [-np.sin(b) * np.cos(a), np.cos(b), -np.sin(b) * np.sin(a)],
                  [-np.sin(a), 0, np.cos(a)]])
    T = np.linalg.inv(T)

    return T
def body_to_inertial(theta, gamma, psi):
    theta = np.radians(theta)
    gamma = np.radians(gamma)
    psi = np.radians(psi)
    T = np.array([[np.cos(theta) * np.cos(psi), np.cos(theta) * np.sin(psi), -np.sin(theta)],
                  [np.sin(gamma) * np.sin(theta) * np.cos(psi) - np.cos(gamma) * np.sin(psi), np.sin(gamma) * np.sin(theta) * np.sin(psi) + np.cos(gamma) * np.cos(psi), np.sin(gamma) * np.sin(theta)],
                  [np.cos(gamma) * np.sin(theta) * np.cos(gamma) + np.sin(gamma) * np.sin(psi), np.cos(gamma) * np.sin(theta) * np.sin(psi) - np.sin(gamma) * np.cos(psi), np.cos(gamma) * np.cos(theta)]])

    T = np.linalg.inv(T)

    return T

    # downwash and instalationa angle are assumed constant
    # aoah = (aoa - aoa0)(1 - deda) + (aoa0 + ih)
    # aoa = (aoal + aoar) /2
    # aoah = aoa - downwash + instalation

#take as input the angle of attack, the speed, the density
#so next thing to do is to get the aerodynamic forces cuz those we cant really coontrol


if __name__ == "__main__":
    R = definegeometry()
    aoa = [7, 7, 9, 7] #left wing, right wing, tail, pitch #first 3 include velocity and rotation aoa, last is pitch for trasformation matrix
    v = [110, 110, 110] #left wing, right wing, tail
    vcruise = [110, 0, 0] #x, y ,z component of velocity
    beta = 0    # mass in kg
    wl, wr, h = aerodynamic(aoa, v, beta)
    pitch = aoa[3]
    W0 = m * 3.71
    # Rwingleft, Rwingright, Rstabilizer, Raileronleft, Raileronright, Rrudder, Rthrustright, Rthrustleft, Cg
    W = np.array([-W0 * np.sin(np.radians(pitch)), 0, W0 * np.cos(np.radians(pitch))])
    F = np.array([wl, wr, h, np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, 0]), W])
    M = moments(R, F)
    print(M)
    a, ang = accelerations(F, M, I)
    print(a)
    Tl = 0
    Tr = 0
    acc = []
    tllist = []
    trlist = []
    times = []
    velo = []
    for i, time in enumerate(np.arange(0, 1000, 0.1)):
        wl, wr, h = aerodynamic(aoa, v, beta)
        pitch = aoa[3]
        W0 = m * 3.71
        W = np.array([-W0 * np.sin(np.radians(pitch)), 0, W0 * np.cos(np.radians(pitch))])
        F = np.array([wl, wr, h, np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([Tl, 0, 0]), np.array([Tr, 0, 0]),
                      np.array([0, 0, 0]), W])
        M = moments(R, F)
        a, ang = accelerations(F, M, I)
        Tl, Tr = thrust(v, aoa, a, Tl, Tr)
        v[0] = v[0] + a[0] * 0.1
        velo.append(v[0])
        acc.append(a[0])
        tllist.append(Tl)
        trlist.append(Tr)
        times.append(time)

    plt.plot(times, tllist, color='b')
    plt.plot(times, trlist, color='r')
    plt.show()
    plt.plot(times, acc, color='g')
    plt.show()
    plt.plot(times, velo, color='b')
    plt.show()