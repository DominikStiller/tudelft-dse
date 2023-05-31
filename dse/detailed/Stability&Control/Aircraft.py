import numpy as np


def define_geometry():
    """
    Returns: the position of each force with respect to the aircraft body reference frame
    """

    Rwingleft = np.array([1, -10, -1])
    Rwingright = np.array([1, 10, -1])
    Rstabilizer = np.array([-10, 0, -1])
    Raileronleft = np.array([0, -16, -1])
    Raileronright = np.array([0, 16, -1])
    Relevator = np.array([-10, 0, -2])
    Rrudder = np.array([-10, 0, -3])
    Cg = np.array([0, 0, 0])
    Rthrustleft = np.array([2, -17, -1])
    Rthrustright = np.array([2, 17, -1])
    R = np.array([Rwingleft, Rwingright, Rstabilizer, Raileronleft, Raileronright, Relevator, Rrudder, Rthrustleft, Rthrustright, Cg])
    return R


def define_areas():
    # plane design parameters
    downwash = 3  # in deg
    instalation = 0  # in deg
    aoastall = 15  # in deg
    area = [56, 56, 15]  # left wing, right wing, stabilizer
    imain = 4  # angle of incidence main
    itail = 4  # angle of incidence tail
    m = 3000
    I = [81250, 19750, 91500]  # mass moments of inertia ixx, iyy, izz
    W0 = m * 3.71
    return area, imain, itail, m, I, W0


def get_coefficients(aoal, aoar, aoah):
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

    aoalistw = np.radians(np.arange(-5, 15, 0.1))
    aoalisth = np.radians(np.arange(-5, 15, 0.1))
    cllistw = aoalistw * 2 * np.pi   # will later be a csv file
    cllisth = aoalisth * 2 * np.pi   # will later be a csv file
    cldivcdw = np.ones(len(aoalistw)) * 20   # will later be a csv file
    cldivcdh = np.ones(len(aoalisth)) * 20   # will later be a csv file
    clwl = np.interp(aoal, aoalistw, cllistw)
    clwr = np.interp(aoar, aoalistw, cllistw)
    clh = np.interp(aoah, aoalisth, cllisth)
    cdwl = clwl / np.interp(aoal, aoalistw, cldivcdw)
    cdwr = clwr / np.interp(aoar, aoalistw, cldivcdw)
    cdh = clh / np.interp(aoah, aoalisth, cldivcdh)
    # add the fuselage contribution later here and add s fuselage
    return clwl, clwr, clh, cdwl, cdwr, cdh

