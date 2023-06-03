import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("TkAgg")

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
    Rthrustleft = np.array([0, -17, -3])
    Rthrustright = np.array([0, 17, -3])
    R = np.array([Rwingleft, Rwingright, Rstabilizer, Raileronleft, Raileronright, Relevator, Rrudder, Rthrustleft, Rthrustright, Cg])
    return R

def define_areas():
    # plane design parameters
    downwash = 3  # in deg
    instalation = 0  # in deg
    aoastall = 15  # in deg
    area = [56, 56, 15]  # left wing, right wing, stabilizer
    imain = np.radians(4)  # angle of incidence main
    itail = np.radians(4)  # angle of incidence tail
    m = 2700
    I = np.array([[81250, 0, 0], [0, 19750, 0], [0, 0, 91500]])  # mass moments of inertia ixx, iyy, izz
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
    if np.radians(-5) <= aoal <= np.radians(15):
        clwl = np.interp(aoal, aoalistw, cllistw)
        cdwl = clwl / np.interp(aoal, aoalistw, cldivcdw)
    else:
        clwl = 0
        cdwl = 1.28
    if np.radians(-5) <= aoar <= np.radians(15):
        clwr = np.interp(aoar, aoalistw, cllistw)
        cdwr = clwr / np.interp(aoar, aoalistw, cldivcdw)
    else:
        clwr = 0
        cdwr = 1.28
    if np.radians(-5) <= aoal <= np.radians(15):
        clh = np.interp(aoah, aoalistw, cllisth)
        cdh = clh / np.interp(aoah, aoalisth, cldivcdh)
    else:
        clh = 0
        cdh = 1.28
    # add the fuselage contribution later and add s fuselage
    return clwl, clwr, clh, cdwl, cdwr, cdh

class System:
    def __init__(self):
        self.geometry = define_geometry()
        self.area, self.imain, self.itail, self.m, self.I, self.W0 = define_areas()
        self.rho = 0.01

        # previous state
        self.euler_prev = np.array([0, 0, 0.])  # the euler angles: roll(at X), pitch(at Y), yaw(at Z)
        self.velocity_linear_prev = np.array([0, 0, 0.])
        self.velocity_angular_prev = np.array([0, 0, 0.])

        # current state
        self.euler = np.copy(self.euler_prev)
        self.velocity_linear = np.copy(self.velocity_linear_prev)
        self.velocity_angular = np.copy(self.velocity_angular_prev)

    def get_state(self):
        return self.euler, self.velocity_linear, self.velocity_angular

    def __call__(self, exc_sig, dt):
        # Fal = [0, 0, 0]
        # Far = [0, 0, 0]
        # Fr = [0, 0, 0]
        Tl, Tr = exc_sig
        aoawings = np.arctan2(self.velocity_linear[2] - self.velocity_angular[1], self.velocity_linear[0])
        aoatail = np.arctan2(self.velocity_linear[2] + self.velocity_angular[1]*10, self.velocity_linear[0])
        # aoa = np.arctan2(self.velocity_linear[2], self.velocity_linear[0])  # if you are ambitious, go do something
        # velocitywings = np.array([self.velocity_linear[0], self.velocity_linear[0], self.velocity_linear[0]]) # if you are ambitious, go do something
        velocitywings = np.array([self.velocity_linear[0] - self.velocity_angular[1], self.velocity_linear[0] - self.velocity_angular[1], self.velocity_linear[0] + self.velocity_angular_prev[1] * 10])
        wl, wr, h = self.aero_forces(aoawings, aoatail, velocitywings, self.euler[2], self.euler[1])
        W = np.array([-self.W0 * np.sin(self.euler[1]), 0, self.W0 * np.cos(self.euler[1])])
        F = np.array([wl, wr, h, np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, 0]), Tl,
                      Tr, W])

        M, Fnet = self.moments(self.geometry, F)
        a, anga = self.accelerations(F, M)

        matrix = np.array([[np.cos(self.euler[1]), np.sin(self.euler[0]) * np.sin(self.euler[1]), np.cos(self.euler[0]) * np.sin(self.euler[1])],
                           [0, np.cos(self.euler[0]) * np.cos(self.euler[1]), -np.sin(self.euler[0]) * np.cos(self.euler[1])],
                           [0, np.sin(self.euler[0]), np.cos(self.euler[0])]])
        d_euler = (1/np.cos(self.euler[1])) * np.matmul(matrix, self.velocity_angular)
        self.euler = self.euler_prev + d_euler*dt
        self.velocity_linear = self.velocity_linear_prev + a*dt
        self.velocity_angular = self.velocity_angular_prev + anga*dt
        self.euler_prev = np.copy(self.euler)
        self.velocity_linear_prev = np.copy(self.velocity_linear)
        self.velocity_angular_prev = np.copy(self.velocity_angular)

        return self.euler, self.velocity_linear, self.velocity_angular

    def aero_forces(self, aoaw, aoah, v, beta, pitch):
        '''
            Args:
                aoaw: angle of attack of the left wing, right wing
                aoah: angle of attack horizontal stabilizer
                v: speed of left wing, right wing and horizontal stabilizer
                beta: sideslip angle
            Returns:
                wl: force vector left wing
                wr: force vector right wing
                h: force vector horizontal stabilizer
        '''

        # aoa includes the effective angle of attack due to the pitch and angular velocity
        clwl, clwr, clh, cdwl, cdwr, cdh = get_coefficients(aoaw + self.imain, aoaw + self.imain, aoah + self.itail)
        Lwl = 0.5 * self.rho * self.area[0] * clwl * v[0] ** 2
        Lwr = 0.5 * self.rho * self.area[1] * clwr * v[1] ** 2
        Lh = 0.5 * self.rho * self.area[2] * clh * v[2] ** 2
        Dwl = 0.5 * self.rho * self.area[0] * cdwl * v[0] ** 2
        Dwr = 0.5 * self.rho * self.area[1] * cdwr * v[1] ** 2
        Dh = 0.5 * self.rho * self.area[2] * cdh * v[2] ** 2

        # Dwl = Dwl + 0.5 * self.rho * self.area[0] * 1 * self.velocity_linear[2] ** 2
        # Dwr = Dwr + 0.5 * self.rho * self.area[1] * 1 * self.velocity_linear[2] ** 2
        # Dh = Dh + 0.5 * self.rho * self.area[2] * 1 * self.velocity_linear[2] ** 2

        wl = np.array([-Dwl, 0, -Lwl])
        wr = np.array([-Dwr, 0, -Lwr])
        h = np.array([-Dh, 0, -Lh])

        wl = np.matmul(self.aero_to_body(pitch, beta), wl)
        wr = np.matmul(self.aero_to_body(pitch, beta), wr)
        h = np.matmul(self.aero_to_body(pitch, beta), h)

        return wl, wr, h

    def aero_to_body(self, aoa, b):
        '''
        Args:
            aoa: angle of attack
            b: angle of sideslip
        Returns:
            T: the transformation matrix
        '''
        a = np.radians(aoa)
        b = np.radians(b)
        T = np.array([[np.cos(b) * np.cos(a), np.sin(b), np.cos(b) * np.sin(a)],
                      [-np.sin(b) * np.cos(a), np.cos(b), -np.sin(b) * np.sin(a)],
                      [-np.sin(a), 0, np.cos(a)]])
        T = np.linalg.inv(T)

        return T

    def moments(self, R, F):
        '''
        Args:
            R: list of the distances for left wing, right wing, horizontal stabilizer, aileron left, aileron right, rudder and cg
            F: list of forces for the same
        Returns:
            M: sum of moment around the X axis, Y axis and Z axis
        '''

        M = np.zeros((len(R), 3))
        for i in range(len(R)):
            M[i] = np.cross(R[i], F[i])
        M = np.array([sum(M.T[0]), sum(M.T[1]), sum(M.T[2])])
        Fnet = np.array([sum(F.T[0]), sum(F.T[1]), sum(F.T[2])])
        return M, Fnet

    def accelerations(self, F, M):
        '''
        Args:
             F: force vector
             M: moment vector
        Returns:
            a: acceleration vector
            ang: angular acceleration vector
        '''

        inv = np.linalg.inv(np.copy(self.I))
        temp = np.matmul(self.I, self.velocity_angular)
        rest = M - np.cross(self.velocity_angular, temp)
        dang_dt = np.matmul(inv, rest)
        temp2 = np.cross(self.velocity_angular, self.velocity_linear)
        a = np.array([sum(F.T[0]) / self.m, sum(F.T[1]) / self.m, sum(F.T[2]) / self.m]) - temp2
        return a, dang_dt

    def body_to_inertial(self):
        '''
        Args:
        Returns:
            T: the transformation matrix
        '''

        theta, gamma, psi = self.euler[1], self.euler[0], self.euler[2]

        T = np.array([[np.cos(theta) * np.cos(psi), np.cos(theta) * np.sin(psi), -np.sin(theta)],
                      [np.sin(gamma) * np.sin(theta) * np.cos(psi) - np.cos(gamma) * np.sin(psi),
                       np.sin(gamma) * np.sin(theta) * np.sin(psi) + np.cos(gamma) * np.cos(psi),
                       np.sin(gamma) * np.sin(theta)],
                      [np.cos(gamma) * np.sin(theta) * np.cos(gamma) + np.sin(gamma) * np.sin(psi),
                       np.cos(gamma) * np.sin(theta) * np.sin(psi) - np.sin(gamma) * np.cos(psi),
                       np.cos(gamma) * np.cos(theta)]])

        T = np.linalg.inv(T)
        return T

class PID:
    def __init__(self, Kp=1., Ki=1., Kd=1.):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd

        self.prev_error = 0.
        self.integral = 0.

    def __call__(self, error, dt):

        # Calculate the integral term
        self.integral += error * dt

        # Calculate the derivative term
        derivative = (error - self.prev_error) / dt

        # Update the previous error for the next iteration
        self.prev_error = error

        return self.Kp * error + self.Ki * self.integral + self.Kd * derivative


class Controller:
    def __init__(self, dt=0.01, Kp=1., Ki=1., Kd=1.):
        self.PID_theta = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_theta):
        return self.PID_theta(error_theta, self.dt)

class Controller2:
    def __init__(self, dt=0.01, Kp=1., Ki=1., Kd=1.):
        self.PID_thrust = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_velocity):
        T = self.PID_thrust(error_velocity, self.dt)
        return T/2, T/2

dt = 0.01
system = System()
controller_theta = Controller(dt=dt, Kp=25000., Ki=4000., Kd=500.)
controller_thrust = Controller2(dt=dt, Kp=1000., Ki=100., Kd=0.)


euler_ref = np.array([0, 0, 0.])
velocity_linear_ref = np.array([0, 0, -2.])
euler_in_time = []
velocity_linear_in_time = []

for i in range(int(1e4)):
    euler, velocity_linear, velocity_angular = system.get_state()

    euler_in_time.append(np.copy(euler))
    velocity_linear_in_time.append(np.copy(velocity_linear))

    error_euler = euler_ref - euler
    error_velocity_linear = velocity_linear_ref - velocity_linear
    print(error_velocity_linear)

    Tl, Tr = controller_thrust(error_velocity_linear)
    # if Tlx>430:
    #     Tlx= 430
    #     Trx= 430
    print(Tl, Tr)
    exc_sig = Tl, Tr

    system(exc_sig, dt)

euler_in_time = np.array(euler_in_time)
velocity_linear_in_time = np.array(velocity_linear_in_time)

# plt.plot(velocity_linear_in_time[:, 0])
# plt.ylabel("velocity in the x direction")
# plt.show()
plt.plot(velocity_linear_in_time[:, 2])
plt.ylabel("velocity in the z direction")
plt.show()
# plt.plot(velocity_linear_in_time[:, 1])
# plt.ylabel("velocity in the y direction")
# plt.show()
# plt.plot(euler_in_time[:, 2], color='b')
# plt.plot(euler_in_time[:, 1], color='g')
# plt.plot(euler_in_time[:, 0], color='r')
# plt.legend('ypr')
# plt.show()
