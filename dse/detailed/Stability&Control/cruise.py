import numpy as np
import scipy as sp

# atmospheric conditions
rho = 0.02
mass = 3000
g = 3.71


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
    Raileronleft = np.array([1, -20, -1])
    Raileronright = np.array([1, 20, -1])
    Rrudder = np.array([-10, 0, -3])
    Cg = np.array([0, 0, 0])
    R = np.array([Rwingleft, Rwingright, Rstabilizer, Raileronleft, Raileronright, Rrudder, Cg])
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


def aerodynamic(aoa, v):
    '''

    Args:
        aoa: angle of attack of the left wing, right wing and horizontal stabilizer
        v: speed of left wing, right wing and horizontal stabilizer

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
    Dwl = 0.5 * rho * area[0] * cdwl * v[0] ** 2
    Dwr = 0.5 * rho * area[1] * cdwr * v[1] ** 2
    Dh = 0.5 * rho * area[2] * cdh * v[2] ** 2
    wl = np.array([-Dwl, 0, Lwl])
    wr = np.array([-Dwr, 0, Lwr])
    h = np.array([-Dh, 0, Lh])
    return wl, wr, h

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


    # downwash and instalationa angle are assumed constant
    # aoah = (aoa - aoa0)(1 - deda) + (aoa0 + ih)
    # aoa = (aoal + aoar) /2
    # aoah = aoa - downwash + instalation

#take as input the angle of attack, the speed, the density
#so next thing to do is to get the aerodynamic forces cuz those we cant really coontrol


if __name__ == "__main__":
    R = definegeometry()
    aoa = [10, 10, 12]
    v = [110, 110, 110]
    wl, wr, h = aerodynamic(aoa, v)
    F = np.array([wl, wr, h, np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, 3000*3.71])])
    moments(R, F)
