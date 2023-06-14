import numpy as np


def define_geometry():
    """
    Returns: the position of each force with respect to the aircraft body reference frame
    """

    # Cessna 172 values
    Rwingleft = np.array([0.259, -4, -1])
    Rwingright = np.array([0.259, 4, -1])
    Rstabilizer = np.array([-3.771, 0, -2])
    Raileronleft = np.array([0, -5, -1])
    Raileronright = np.array([0, 5, -1])
    Relevator = np.array([-3.771, 0, -1])
    Rrudder = np.array([-3.771, 0, -2])
    Cg = np.array([0, 0, 0])
    Rthrustleft = np.array([3, 0, 0])
    Rthrustright = np.array([0, 0, -0])

    # Rwingleft = np.array([0.356, -10, -1])
    # Rwingright = np.array([0.356, 10, -1])
    # Rstabilizer = np.array([-12.644, 0, -1])
    # Raileronleft = np.array([0.356, -20, -1])
    # Raileronright = np.array([0.356, 20, -1])
    # Relevator = np.array([-12.644, 0, -1])
    # Rrudder = np.array([-12.644, 0, -3])
    # Cg = np.array([0, 0, 0])
    # Rthrustleft = np.array([1.856, -22, -1])
    # Rthrustright = np.array([1.856, 22, -1])
    R = np.array([Rwingleft, Rwingright, Rstabilizer, Raileronleft, Raileronright, Relevator, Rrudder, Rthrustleft, Rthrustright, Cg])
    return R


def define_areas():
    # # plane design parameters
    # downwash = 3  # in deg
    # instalation = 0  # in deg
    # aoastall = 15  # in deg
    # area = [66.5, 66.5, 16.9]  # left wing, right wing, stabilizer
    # imain = np.radians(1)  # angle of incidence main 10
    # itail = np.radians(2)  # angle of incidence tail
    # m = 3000
    # I = np.array([[98156, 0, 0], [0, 5645, 0], [0, 0, 98310]])
    # W0 = m * 3.71

    # Cessna parameters
    area = [8.1, 8.1, 16.2*0.202]  # left wing, right wing, stabilizer
    imain = np.radians(0)
    itail = np.radians(0)
    I = np.array([[1285, 0, 0], [0, 1825, 0], [0, 0, 2667]])
    m = 1043
    W0 = m * 9.81
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

    # aoalistw = np.radians(np.arange(-180, 181, 5))
    # aoalisth = np.radians(np.arange(-180, 181, 5))

    # cllistw = np.array(
    #     [0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
    #      -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.6, -0.8, -0.1, 0.2, 0.6,
    #      1.2, 1.9, 0.8, 0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
    #      -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5, -0.1])
    # cllisth = np.array(
    #     [-0.1, 0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
    #      -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.6, -0.8, -0.1, 0.2, 0.6,
    #      1.2, 1.9, 0.8, 0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
    #      -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5])

    # cdwr = 0.65 - 0.64 * np.cos(2 * aoar)
    # cdwl = 0.65 - 0.64 * np.cos(2 * aoal)
    # cdh = 0.65 - 0.64 * np.cos(2 * aoah)

    # clwl = np.interp(aoal, aoalistw, cllistw)
    # clwr = np.interp(aoar, aoalistw, cllistw)
    # clh = np.interp(aoah, aoalisth, cllisth)

    # Our plane parameters
    # if aoal > np.radians(15) or aoal < np.radians(-5):
    #     aoalistw = np.radians(np.arange(-180, 181, 5))
    #     cllistw = np.array(
    #         [0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
    #          -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.6, -0.2, 0.4, 1.2, 1.75,
    #          2.2, 2.3, 1.5, 0.8, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
    #          -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5, -0.1])
    #     clwl = np.interp(aoal, aoalistw, cllistw)
    # else:
    #     aoalistw = np.radians(np.arange(-5, 15.1, 1))
    #     cllistw = np.array([0.4, 0.4, 0.5, 0.9, 1.1, 1.2, 1.33, 1.43, 1.53, 1.64, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.3, 2.35, 2.4, 2.35, 2.3])
    #     clwl = np.interp(aoal, aoalistw, cllistw)
    #
    # if aoar > np.radians(15) or aoar < np.radians(-5):
    #     aoalistw = np.radians(np.arange(-180, 181, 5))
    #     cllistw = np.array(
    #         [0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
    #          -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.6, -0.2, 0.4, 1.2, 1.75,
    #          2.2, 2.3, 1.5, 0.8, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
    #          -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5, -0.1])
    #     clwr = np.interp(aoar, aoalistw, cllistw)
    # else:
    #     aoalistw = np.radians(np.arange(-5, 15.1, 1))
    #     cllistw = np.array([0.4, 0.4, 0.5, 0.9, 1.1, 1.2, 1.33, 1.43, 1.53, 1.64, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.3, 2.35, 2.4, 2.35, 2.3])
    #     clwr = np.interp(aoar, aoalistw, cllistw)
    #
    # if aoah > np.radians(15) or aoah < np.radians(-5):
    #     aoalisth = np.radians(np.arange(-180, 181, 5))
    #     cllisth = np.array(
    #         [0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
    #          -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, -0.7, -1.2, -1, -0.6, 0.0, 0.6,
    #          1, 1.5, 0.8, 0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
    #          -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5, -0.1])
    #     clh = np.interp(aoah, aoalisth, cllisth)
    # else:
    #     aoalisth = np.radians(np.arange(-5, 15.1, 1))
    #     cllisth = np.array([-0.6, -0.5, -0.4, -0.3, -0.15, 0, 0.15, 0.327, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.5])
    #     clh = np.interp(aoah, aoalisth, cllisth)
    #
    #
    # cdwr = 0.65 - 0.64 * np.cos(2 * aoar)
    # cdwl = 0.65 - 0.64 * np.cos(2 * aoal)
    # cdh = 0.65 - 0.64 * np.cos(2 * aoah)


    # Cessna 172 Parameters
    if aoal > np.radians(15) or aoal < np.radians(-5):
        aoalistw = np.radians(np.arange(-180, 181, 5))
        cllistw = np.array(
            [0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
             -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.6, -0.2, -0.4, 1.2, 1.75,
             2.2, 2.3, 1.5, 0.8, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
             -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5, -0.1])
        clwl = np.interp(aoal, aoalistw, cllistw)
    else:
        aoalistw = np.radians(np.arange(-5, 15.1, 1))
        cllistw = np.array([-0.4, -0.3, -0.2, -0.05, 0.1, 0.25, 0.327, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.3, 1.3, 1.3, 1.35, 1.4])
        clwl = np.interp(aoal, aoalistw, cllistw)

    if aoar > np.radians(15) or aoar < np.radians(-5):
        aoalistw = np.radians(np.arange(-180, 181, 5))
        cllistw = np.array(
            [0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
             -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.6, -0.2, 0.4, 1.2, 1.75,
             2.2, 2.3, 1.5, 0.8, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
             -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5, -0.1])
        clwr = np.interp(aoar, aoalistw, cllistw)
    else:
        aoalistw = np.radians(np.arange(-5, 15.1, 1))
        cllistw = np.array([-0.4, -0.3, -0.2, -0.05, 0.1, 0.25, 0.327, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.3, 1.3, 1.3, 1.35, 1.4])
        clwr = np.interp(aoar, aoalistw, cllistw)

    if aoah > np.radians(15) or aoah < np.radians(-5):
        aoalisth = np.radians(np.arange(-180, 181, 5))
        cllisth = np.array(
            [0.15, 0.4, 0.6, 0.65, 0.2, 0.4, 0.7, 0.9, 1.1, 1.2, 1.1, 1.0, 0.95, 0.8, 0.7, 0.55, 0.3, 0.0, -0.3, -0.45,
             -0.65, -0.75, -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, -0.7, -1.2, -1, -0.6, 0.0, 0.6,
             1, 1.5, 0.8, 0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.15, 1.1, 1.0, 0.8, 0.6, 0.4, 0.2, 0, -0.3, -0.45, -0.65, -0.75,
             -0.85, -0.95, -1.05, -1.1, -1.15, -1.2, -1.0, -0.8, -0.6, -0.3, 0, -0.3, -0.5, -0.1])
        clh = np.interp(aoah, aoalisth, cllisth)
    else:
        aoalisth = np.radians(np.arange(-5, 15.1, 1))
        cllisth = np.array([-0.6, -0.5, -0.4, -0.3, -0.15, 0, 0.156, 0.327, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.5])
        clh = np.interp(aoah, aoalisth, cllisth)


    cdwr = 0.65 - 0.64 * np.cos(2 * aoar)
    cdwl = 0.65 - 0.64 * np.cos(2 * aoal)
    cdh = 0.65 - 0.64 * np.cos(2 * aoah)

    return clwl, clwr, clh, cdwl, cdwr, cdh

