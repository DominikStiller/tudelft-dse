import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("MacOSX")

def define_geometry(ithrust):
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
    # print(f"{np.degrees(ithrust)}---------")
    Rthrustleft = np.array([0+3*np.cos(ithrust[-1]), -17, -3*np.sin(ithrust[-1])])
    Rthrustright = np.array([0+3*np.cos(ithrust[-1]), 17, -3*np.sin(ithrust[-1])])
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
    max_thrust = -11000
    return area, imain, itail, m, I, W0, max_thrust

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

    cdwr = 0.65 - 0.64 * np.cos(2 * aoar)
    cdwl = 0.65 - 0.64 * np.cos(2 * aoal)
    cdh = 0.65 - 0.64 * np.cos(2 * aoah)

    clwl = np.interp(aoal, aoalistw, cllistw)
    clwr = np.interp(aoar, aoalistw, cllistw)
    clh = np.interp(aoah, aoalisth, cllisth)

    # add the fuselage contribution later and add s fuselage
    return clwl, clwr, clh, cdwr, cdwl, cdh

class System:
    def __init__(self):
        self.ithrust_prev = np.array([np.radians(90)])
        self.geometry = define_geometry(self.ithrust_prev)
        self.area, self.imain, self.itail, self.m, self.I, self.W0, self.max_thrust = define_areas()
        self.rho = 0.01

        # previous state
        self.euler_prev = np.array([0, 0, 0.])  # the euler angles: roll(at X), pitch(at Y), yaw(at Z)
        self.velocity_linear_prev = np.array([0, 0, 0])
        self.velocity_angular_prev = np.array([0, 0, 0.])
        self.position_prev = np.array([0, 0, 0])

        # current state
        self.euler = np.copy(self.euler_prev)
        self.velocity_linear = np.copy(self.velocity_linear_prev)
        self.velocity_angular = np.copy(self.velocity_angular_prev)
        self.position = np.copy(self.position_prev)
        self.ithrust = np.copy(self.ithrust_prev)

    def get_state(self):
        return self.euler, self.velocity_linear, self.velocity_angular, self.position, self.ithrust

    def __call__(self, exc_sig, dt):
        # Fal = [0, 0, 0]
        # Far = [0, 0, 0]
        # Fr = [0, 0, 0]
        Tl, Tr, Fele = exc_sig

        ithrust_prev = self.ithrust_prev
        ithrust_current = np.array([np.arctan2(2*Tl[2], 2*Tl[0])])
        self.ithrust = ithrust_current
        # print(f"thrust ang prev: {ithrust_prev}, thrust ang current: {ithrust_current}")
        # if ithrust_prev - ithrust_current >= np.radians(3) * dt:
        #     self.ithrust = ithrust_prev - np.radians(3) * dt
        # else:
        #     self.ithrust = ithrust_current
        #
        # if ithrust_prev <= 0:
        #     self.ithrust = ithrust_prev + np.radians(3) * dt
        # if np.sqrt((2 * Tl[0]) ** 2 + (2 * Tl[1]) ** 2 + (2 * Tl[2]) ** 2) > abs(system.max_thrust):
        #     Tl[2] = (W[2] + wl[2] + wr[2] + h[2]) / 2
        #     Tr[2] = (W[2] + wl[2] + wr[2] + h[2]) / 2
        #     Tl[0] = np.sqrt(system.max_thrust ** 2 - 2*Tl[2] ** 2) / 2
        #     Tr[0] = np.sqrt(system.max_thrust ** 2 - 2*Tr[2] ** 2) / 2
        # else:
        #     Tl[0] = (Tl[2] / np.tan(self.ithrust))[0]
        #     Tr[0] = (Tr[2] / np.tan(self.ithrust))[0]

        self.geometry = define_geometry(self.ithrust)
        aoawings = np.arctan2(self.velocity_linear[2] - self.velocity_angular[1], self.velocity_linear[0])
        print(aoawings, "AoA wing")
        aoatail = np.arctan2(self.velocity_linear[2] + self.velocity_angular[1]*10, self.velocity_linear[0])
        velocitywings = np.array([self.velocity_linear[0] - self.velocity_angular[1], self.velocity_linear[0] - self.velocity_angular[1], self.velocity_linear[0] + self.velocity_angular_prev[1] * 10])
        wl, wr, h = self.aero_forces(aoawings, aoatail, velocitywings, self.euler[2], self.euler[1])
        W = np.array([-self.W0 * np.sin(self.euler[1]), 0, self.W0 * np.cos(self.euler[1])])

        # print(f"Lift force of wing = {wl+wr}")
        # print(f"Lift force of horizontal stabilizer = {h}")
        # print(f"Elevator force = {Fele}")
        # print(f"Thrust force = {Tl+Tr}")
        # print(f"{np.sqrt((2*Tl[0])**2 + (2*Tl[1])**2 + (2*Tl[2])**2)}---------")

        self.F = np.array([wl, wr, h, np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, Fele]), np.array([0, 0, 0]), Tl,
                      Tr, W])

        M, Fnet = self.moments(self.geometry, self.F)
        a, anga = self.accelerations(self.F, M)
        anga = 0
        matrix = np.array([[np.cos(self.euler[1]), np.sin(self.euler[0]) * np.sin(self.euler[1]), np.cos(self.euler[0]) * np.sin(self.euler[1])],
                           [0, np.cos(self.euler[0]) * np.cos(self.euler[1]), -np.sin(self.euler[0]) * np.cos(self.euler[1])],
                           [0, np.sin(self.euler[0]), np.cos(self.euler[0])]])
        d_euler = (1/np.cos(self.euler[1])) * np.matmul(matrix, self.velocity_angular)
        self.euler = self.euler_prev + d_euler*dt
        self.velocity_linear = self.velocity_linear_prev + a*dt
        self.velocity_angular = self.velocity_angular_prev + anga*dt
        self.position = self.position_prev + self.velocity_linear*dt

        self.euler_prev = np.copy(self.euler)
        self.velocity_linear_prev = np.copy(self.velocity_linear)
        self.velocity_angular_prev = np.copy(self.velocity_angular)
        self.position_prev = np.copy(self.position)
        self.ithrust_prev = np.copy(self.ithrust)

        return self.euler, self.velocity_linear, self.velocity_angular, self.position, self.ithrust

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

        Dwlz = 0.5 * self.rho * self.area[0] * cdwlz * self.velocity_linear[2] ** 2 * -np.sign(self.velocity_linear[2])
        Dwrz = 0.5 * self.rho * self.area[1] * cdwrz * self.velocity_linear[2] ** 2 * -np.sign(self.velocity_linear[2])
        Dhz = 0.5 * self.rho * self.area[2] * cdhz * self.velocity_linear[2] ** 2 * -np.sign(self.velocity_linear[2])

        wl = np.array([Dwlx, 0, -Lwl + Dwlz])
        wr = np.array([Dwrx, 0, -Lwr + Dwrz])
        h = np.array([Dhx, 0, -Lh + Dhz])

        wl = np.matmul(self.aero_to_body(aoaw, beta), wl)
        wr = np.matmul(self.aero_to_body(aoaw, beta), wr)
        h = np.matmul(self.aero_to_body(aoah, beta), h)

        return wl, wr, h

    def aero_to_body(self, aoa, b):
        '''
        Args:
            aoa: angle of attack
            b: angle of sideslip
        Returns:
            T: the transformation matrix
        '''
        a = aoa
        b = b
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
    def __init__(self, dt=0.01, Kp=0., Ki=0., Kd=0.):
        self.PID_zthrust = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_velocity):
        T = self.PID_thrust(error_velocity, self.dt)
        # if T < system.max_thrust:
        #     T = system.max_thrust
        #     return T / 2, T / 2
        # elif T > 0:
        #     T = 0
        #     return T / 2, T / 2
        # else:
        return T / 2, T / 2

class Controller2:
    def __init__(self, dt=0.01, Kp=0., Ki=0., Kd=0.):
        self.PID_thrust_pitch = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, pitch):
        T = self.PID_thrust_pitch(pitch, self.dt)
        return T / 2, T / 2

class Controller3:
    def __init__(self, dt=0.01, Kp=0., Ki=0., Kd=0.):
        self.PID_xthrust = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, xthrust):
        ithrust = self.PID_xthrust(xthrust, self.dt)
        return ithrust / 2, ithrust / 2

class Controller4:
    def __init__(self, dt=0.01, Kp=1., Ki=1., Kd=1.):
        self.PID_theta = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_theta):
        ele = self.PID_theta(error_theta, self.dt)
        return ele

# def TransitionSimple(V, tilt, cl, S, clWING):
#     #SOURCE: PRouty
#     R = rotorParameters["Radius"]
#     rho = const["airDensity"]
#     b= rotorParameters["N_blades"]
#     Vmax = const["soundSpeed"]*0.92
#     powmax = rotorParameters['maxEnginePower']
#     Vtang = V*np.cos(tilt)
#     omega = (Vmax-Vtang)/R
#     mu = V/(omega*R) * np.cos(tilt)
#     c = R/20
#     func = lambda r, phi: (omega * R * (r / R + mu * np.sin(phi)))**2 * rho/2 * cl * c
#     T = b/(2*np.pi) * integrate.dblquad(func, 0, 2*np.pi, 0, R)[0]
#     if T*np.sin(tilt)*V > powmax:
#         T = powmax/(V*np.sin(tilt))
#     L = 1/2 * rho * (V+44*np.sin(tilt))**2 * S * clWING
#     Tvert = T*np.cos(tilt)*4+L
#
#     return T*4, Tvert

class Controller5:
    def __init__(self, dt=0.01, Kp=0., Ki=0., Kd=0.):
        self.PID_thrust = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_mag_thrust, T_current, error_velocity_linear):
        T_change = self.PID_thrust(error_mag_thrust, self.dt)
        if abs(T_change) > abs(system.max_thrust):
            T_change = system.max_thrust
        # print(T_change, "tchangetchangetchange")

        # T_target = abs(np.linalg.norm(T_current) + T_change)
        T = np.array([0, 0, 0])
        T[2] = -(system.F[9][2] + system.F[0][2] + system.F[1][2] + system.F[2][2])
        print(system.F[9][2], "Weight")
        print(system.F[0][2], "Left wing")
        print(system.F[1][2], "Right wing")
        print(system.F[2][2], "Horizontal stabilizer")
        controller_xthrust = Controller3(dt=dt, Kp=10., Ki=0., Kd=5.)
        Tlx, Trx = controller_xthrust(error_velocity_linear[0])
        if abs(error_velocity_linear[0]) < 20:
            T[0] = Tlx + Trx
        else:
            T[0] = np.sqrt(T_change ** 2 - T[2] ** 2)
        # if np.linalg.norm(T_target) > abs(system.max_thrust):
        #     T[2] = (W[2] + wl[2] + wr[2] + h[2])
        #     T[0] = np.sqrt(system.max_thrust ** 2 - 2 * Tl[2] ** 2)
        # else:
        #     T[0] = (T[2] / np.tan(ithrust))[0]
        return T

dt = 0.01
system = System()
controller_thrust = Controller(dt=dt, Kp=2000., Ki=0., Kd=300.)
# controller_pitch = Controller2(dt=dt, Kp=0., Ki=0., Kd=0.)
controller_xthrust = Controller3(dt=dt, Kp=6000., Ki=2500., Kd=3000.)
controller_elevator = Controller4(dt=dt, Kp=12000., Ki=600., Kd=4000.)
controller_control_thrust = Controller5(dt=dt, Kp=1000., Ki=0., Kd=2500.)

euler_ref = np.array([0., 0., 0.])
velocity_linear_ref = np.array([400/3.6, 0., 0.])
ithrust_ref = np.array([0])
control_thrust_ref = system.max_thrust

euler_in_time = []
velocity_linear_in_time = []
position_in_time = []
thrust_in_time = []
ithrust_in_time = []

# Initial state call
Tl, Tr = np.array([0, 0, -5500]), np.array([0, 0, -5500])
Fele = 0
exc_sig = Tl, Tr, Fele
system(exc_sig, dt)

for i in range(int(1e4)):
    euler, velocity_linear, velocity_angular, position, ithrust = system.get_state()

    euler_in_time.append(np.copy(euler))
    print(euler[1], "pitch")
    velocity_linear_in_time.append(np.copy(velocity_linear))
    position_in_time.append(np.copy(position))
    ithrust_in_time.append(np.copy(ithrust))

    error_euler = euler_ref - euler
    error_velocity_linear = velocity_linear_ref - velocity_linear
    error_ithrust = ithrust_ref - ithrust
    # print(error_velocity_linear)

    Tlz, Trz = controller_zthrust(error_velocity_linear[2])
    Tlx, Trx = controller_xthrust(error_velocity_linear[0])

    Fele = controller_elevator(error_euler[1])
    # print(Fele, "elevator force")

    Tl = np.array([Tlx, 0, Tlz])
    Tr = np.array([Trx, 0, Trz])
    T =  Tl + Tr
    # print(T, "t1t1t1")

    error_mag_thrust = abs(control_thrust_ref) - np.linalg.norm(T)


    T = controller_control_thrust(error_mag_thrust, T, error_velocity_linear)

    thrust_in_time.append(np.copy(T))

    print(T, velocity_linear, "t2+vt2+vt2+v")
    exc_sig = T / 2, T / 2, Fele

    system(exc_sig, dt)

euler_in_time = np.array(euler_in_time)
velocity_linear_in_time = np.array(velocity_linear_in_time)
position_in_time = np.array(position_in_time)
thrust_in_time = np.array(thrust_in_time)

plt.plot(velocity_linear_in_time[:, 0])
plt.ylabel("velocity in the x direction")
plt.show()
# plt.plot(velocity_linear_in_time[:, 1])
# plt.ylabel("velocity in the y direction")
# plt.show()
plt.plot(velocity_linear_in_time[:, 2])
plt.ylabel("velocity in the z direction")
plt.show()

plt.plot(position_in_time[:, 2], label="z-position")
plt.plot(position_in_time[:, 0], label="x-position")
plt.ylabel("position")
plt.legend()
plt.show()

plt.plot(thrust_in_time[:,2], label="z-thrust")
plt.plot(thrust_in_time[:,0], label="x-thrust")
plt.legend()
plt.show()

plt.plot(np.degrees(ithrust_in_time[:]))
plt.ylabel("thrust angle [deg]")
plt.show()

plt.plot(euler_in_time[:, 2], color='b')
plt.plot(euler_in_time[:, 1], color='g')
plt.plot(euler_in_time[:, 0], color='r')
plt.legend('ypr')
plt.show()