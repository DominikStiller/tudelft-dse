import numpy as np


def define_geometry():
    """
    Returns: the position of each force with respect to the aircraft body reference frame
    """

    Rwingleft = np.array([1, -10, -1])
    Rwingright = np.array([1, 10, -1])
    Rstabilizer = np.array([-10, 0, -2])
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
    imain = np.radians(10)  # angle of incidence main
    itail = np.radians(10)  # angle of incidence tail
    m = 3000
    #I = np.array([[81250, 0, 0], [0, 19750, 0], [0, 0, 91500]])  # mass moments of inertia ixx, iyy, izz
    I = np.array([[98156, 0, 0], [0, 5645, 0], [0, 0, 98310]])
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

    # aoalistw = np.radians(np.arange(-5, 15, 0.1))
    # aoalisth = np.radians(np.arange(-5, 15, 0.1))
    # cllistw = aoalistw * 2 * np.pi   # will later be a csv file
    # cllisth = aoalisth * 2 * np.pi   # will later be a csv file
    # cldivcdw = np.ones(len(aoalistw)) * 20   # will later be a csv file
    # cldivcdh = np.ones(len(aoalisth)) * 20   # will later be a csv file
    # if np.radians(-5) <= aoal <= np.radians(15):
    #     clwl = np.interp(aoal, aoalistw, cllistw)
    #     cdwl = clwl / np.interp(aoal, aoalistw, cldivcdw)
    # else:
    #     clwl = 0
    #     cdwl = 1
    # if np.radians(-5) <= aoar <= np.radians(15):
    #     clwr = np.interp(aoar, aoalistw, cllistw)
    #     cdwr = clwr / np.interp(aoar, aoalistw, cldivcdw)
    # else:
    #     clwr = 0
    #     cdwr = 1
    # if np.radians(-5) <= aoal <= np.radians(15):
    #     clh = np.interp(aoah, aoalistw, cllisth)
    #     cdh = clh / np.interp(aoah, aoalisth, cldivcdh)
    # else:
    #     clh = 0
    #     cdh = 1

    aoalistw = np.radians(np.arange(-180, 181, 5))
    aoalisth = np.radians(np.arange(-180, 181, 5))
    cllistw = np.array(
        [0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
         -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.6, -0.8, -0.1, 0.2, 0.6,
         1.2, 1.9, 0.8, 0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
         -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5, -0.1])
    cllisth = np.array(
        [-0.1, 0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
         -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.6, -0.8, -0.1, 0.2, 0.6,
         1.2, 1.9, 0.8, 0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
         -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5])

    cdwr = 0.65 - 0.64 * np.cos(2 * np.radians(aoar))
    cdwl = 0.65 - 0.64 * np.cos(2 * np.radians(aoal))
    cdh = 0.64 - 0.64 * np.cos(2 * np.radians(aoah))

    clwl = np.interp(aoal, aoalistw, cllistw)
    clwr = np.interp(aoar, aoalistw, cllistw)
    clh = np.interp(aoah, aoalisth, cllisth)
    # add the fuselage contribution later and add s fuselage
    return clwl, clwr, clh, cdwl, cdwr, cdh

