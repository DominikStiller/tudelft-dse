import matplotlib.pyplot as plt
import numpy as np

from dse.plotting import format_plot


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
    Rthrustleft = np.array([0 + 3 * np.cos(ithrust[-1]), -17, -3 * np.sin(ithrust[-1])])
    Rthrustright = np.array([0 + 3 * np.cos(ithrust[-1]), 17, -3 * np.sin(ithrust[-1])])
    R = np.array(
        [
            Rwingleft,
            Rwingright,
            Rstabilizer,
            Raileronleft,
            Raileronright,
            Relevator,
            Rrudder,
            Rthrustleft,
            Rthrustright,
            Cg,
        ]
    )
    return R


def define_areas():
    area = [66.75, 66.75, 16.9, 31]  # left wing, right wing, stabilizer, fuselage
    imain = np.radians(1)  # angle of incidence main
    itail = np.radians(1)  # angle of incidence tail
    m = 3000
    I = np.array([[98156, 0, 0], [0, 5645, 0], [0, 0, 98310]])
    W0 = m * 3.71
    max_thrust = -12284
    return area, imain, itail, m, I, W0, max_thrust


def get_coefficients(aoal, aoar, aoah):
    """
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
    """
    if aoal > np.radians(15) or aoal < np.radians(-5):
        aoalistw = np.radians(np.arange(-180, 181, 5))
        cllistw = np.array(
            [
                0.15,
                0.4,
                0.6,
                0.65,
                0.2,
                0.4,
                0.7,
                0.9,
                1.1,
                1.2,
                1.1,
                1.0,
                0.95,
                0.8,
                0.7,
                0.55,
                0.3,
                0.0,
                -0.3,
                -0.45,
                -0.65,
                -0.75,
                -0.85,
                -0.95,
                -1.05,
                -1.1,
                -1.15,
                -1.2,
                -1.0,
                -0.8,
                -0.6,
                -0.3,
                0,
                -0.6,
                -0.2,
                0.4,
                1.2,
                1.75,
                2.2,
                2.3,
                1.5,
                0.8,
                0.6,
                0.8,
                0.9,
                1,
                1.1,
                1.15,
                1.1,
                1.0,
                0.8,
                0.6,
                0.4,
                0.2,
                0,
                -0.3,
                -0.45,
                -0.65,
                -0.75,
                -0.85,
                -0.95,
                -1.05,
                -1.1,
                -1.15,
                -1.2,
                -1.0,
                -0.8,
                -0.6,
                -0.3,
                0,
                -0.3,
                -0.5,
                -0.1,
            ]
        )
        clwl = np.interp(aoal, aoalistw, cllistw)
    else:
        aoalistw = np.radians(np.arange(-5, 15.1, 1))
        cllistw = np.array(
            [
                0.4,
                0.4,
                0.5,
                0.9,
                1.1,
                1.2,
                1.33,
                1.43,
                1.53,
                1.64,
                1.75,
                1.85,
                1.95,
                2.05,
                2.15,
                2.25,
                2.3,
                2.35,
                2.4,
                2.35,
                2.3,
            ]
        )
        clwl = np.interp(aoal, aoalistw, cllistw)

    if aoar > np.radians(15) or aoar < np.radians(-5):
        aoalistw = np.radians(np.arange(-180, 181, 5))
        cllistw = np.array(
            [
                0.15,
                0.4,
                0.6,
                0.65,
                0.2,
                0.4,
                0.7,
                0.9,
                1.1,
                1.2,
                1.1,
                1.0,
                0.95,
                0.8,
                0.7,
                0.55,
                0.3,
                0.0,
                -0.3,
                -0.45,
                -0.65,
                -0.75,
                -0.85,
                -0.95,
                -1.05,
                -1.1,
                -1.15,
                -1.2,
                -1.0,
                -0.8,
                -0.6,
                -0.3,
                0,
                -0.6,
                -0.2,
                0.4,
                1.2,
                1.75,
                2.2,
                2.3,
                1.5,
                0.8,
                0.6,
                0.8,
                0.9,
                1,
                1.1,
                1.15,
                1.1,
                1.0,
                0.8,
                0.6,
                0.4,
                0.2,
                0,
                -0.3,
                -0.45,
                -0.65,
                -0.75,
                -0.85,
                -0.95,
                -1.05,
                -1.1,
                -1.15,
                -1.2,
                -1.0,
                -0.8,
                -0.6,
                -0.3,
                0,
                -0.3,
                -0.5,
                -0.1,
            ]
        )
        clwr = np.interp(aoar, aoalistw, cllistw)
    else:
        aoalistw = np.radians(np.arange(-5, 15.1, 1))
        cllistw = np.array(
            [
                0.4,
                0.4,
                0.5,
                0.9,
                1.1,
                1.2,
                1.33,
                1.43,
                1.53,
                1.64,
                1.75,
                1.85,
                1.95,
                2.05,
                2.15,
                2.25,
                2.3,
                2.35,
                2.4,
                2.35,
                2.3,
            ]
        )
        clwr = np.interp(aoar, aoalistw, cllistw)

    if aoah > np.radians(15) or aoah < np.radians(-5):
        aoalisth = np.radians(np.arange(-180, 181, 5))
        cllisth = np.array(
            [
                0.15,
                0.4,
                0.6,
                0.65,
                0.2,
                0.4,
                0.7,
                0.9,
                1.1,
                1.2,
                1.1,
                1.0,
                0.95,
                0.8,
                0.7,
                0.55,
                0.3,
                0.0,
                -0.3,
                -0.45,
                -0.65,
                -0.75,
                -0.85,
                -0.95,
                -1.05,
                -1.1,
                -1.15,
                -1.2,
                -1.0,
                -0.8,
                -0.6,
                -0.3,
                -0.7,
                -1.2,
                -1,
                -0.6,
                0.0,
                0.6,
                1,
                1.5,
                0.8,
                0.4,
                0.6,
                0.8,
                0.9,
                1,
                1.1,
                1.15,
                1.1,
                1.0,
                0.8,
                0.6,
                0.4,
                0.2,
                0,
                -0.3,
                -0.45,
                -0.65,
                -0.75,
                -0.85,
                -0.95,
                -1.05,
                -1.1,
                -1.15,
                -1.2,
                -1.0,
                -0.8,
                -0.6,
                -0.3,
                0,
                -0.3,
                -0.5,
                -0.1,
            ]
        )
        clh = np.interp(aoah, aoalisth, cllisth)
    else:
        aoalisth = np.radians(np.arange(-5, 15.1, 1))
        cllisth = np.array(
            [
                -0.6,
                -0.5,
                -0.4,
                -0.3,
                -0.15,
                0,
                0.1635,
                0.327,
                0.4,
                0.5,
                0.6,
                0.7,
                0.8,
                0.9,
                1.0,
                1.1,
                1.2,
                1.3,
                1.4,
                1.5,
                1.5,
            ]
        )
        clh = np.interp(aoah, aoalisth, cllisth)

    cdwr = 0.6625 - 0.6175 * np.cos(2 * aoar)
    cdwl = 0.6625 - 0.6175 * np.cos(2 * aoal)
    cdh = 0.6625 - 0.6175 * np.cos(2 * aoah)
    cdb = 0.031

    return clwl, clwr, clh, cdwr, cdwl, cdh, cdb

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


class System:
    def __init__(self):
        self.ithrust_prev = np.array([np.radians(90)])
        self.geometry = define_geometry(self.ithrust_prev)
        self.area, self.imain, self.itail, self.m, self.I, self.W0, self.max_thrust = define_areas()
        self.rho = 0.01

        # previous state
        self.euler_prev = np.array(
            [0, 0, 0.0]
        )  # the euler angles: roll(at X), pitch(at Y), yaw(at Z)
        self.velocity_linear_prev = np.array([0, 0, 0])
        self.velocity_angular_prev = np.array([0, 0, 0.0])
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
        ithrust_current = np.array([np.arctan2(2 * Tl[2], 2 * Tl[0])])
        self.ithrust = ithrust_current

        self.geometry = define_geometry(self.ithrust)
        aoawings = np.arctan2(
            self.velocity_linear[2] - self.velocity_angular[1], self.velocity_linear[0]
        )
        aoatail = np.arctan2(
            self.velocity_linear[2] + self.velocity_angular[1] * 13, self.velocity_linear[0]
        )
        velocitywings = np.array(
            [
                self.velocity_linear[0] - self.velocity_angular[1],
                self.velocity_linear[0] - self.velocity_angular[1],
                self.velocity_linear[0] + self.velocity_angular_prev[1] * 10,
            ]
        )
        wl, wr, h, b = self.aero_forces(
            aoawings, aoatail, velocitywings, self.euler[2], self.euler[1]
        )
        W = np.array([-self.W0 * np.sin(self.euler[1]), 0, self.W0 * np.cos(self.euler[1])])
        W = W + b

        self.F = np.array(
            [
                wl,
                wr,
                h,
                np.array([0, 0, 0]),
                np.array([0, 0, 0]),
                np.array([0, 0, Fele]),
                np.array([0, 0, 0]),
                Tl,
                Tr,
                W,
            ]
        )

        M, Fnet = self.moments(self.geometry, self.F)
        a, anga = self.accelerations(self.F, M)
        anga = 0
        matrix = np.array(
            [
                [
                    np.cos(self.euler[1]),
                    np.sin(self.euler[0]) * np.sin(self.euler[1]),
                    np.cos(self.euler[0]) * np.sin(self.euler[1]),
                ],
                [
                    0,
                    np.cos(self.euler[0]) * np.cos(self.euler[1]),
                    -np.sin(self.euler[0]) * np.cos(self.euler[1]),
                ],
                [0, np.sin(self.euler[0]), np.cos(self.euler[0])],
            ]
        )
        d_euler = (1 / np.cos(self.euler[1])) * np.matmul(matrix, self.velocity_angular)
        self.euler = self.euler_prev + d_euler * dt
        self.velocity_linear = self.velocity_linear_prev + a * dt
        self.velocity_angular = self.velocity_angular_prev + anga * dt
        self.position = self.position_prev + self.velocity_linear * dt

        self.euler_prev = np.copy(self.euler)
        self.velocity_linear_prev = np.copy(self.velocity_linear)
        self.velocity_angular_prev = np.copy(self.velocity_angular)
        self.position_prev = np.copy(self.position)
        self.ithrust_prev = np.copy(self.ithrust)

        return self.euler, self.velocity_linear, self.velocity_angular, self.position, self.ithrust

    def aero_forces(self, aoaw, aoah, v, beta, pitch):
        """
        Args:
            aoaw: angle of attack of the left wing, right wing
            aoah: angle of attack horizontal stabilizer
            v: speed of left wing, right wing and horizontal stabilizer
            beta: sideslip angle
        Returns:
            wl: force vector left wing
            wr: force vector right wing
            h: force vector horizontal stabilizer
        """

        # aoa includes the effective angle of attack due to the pitch and angular velocity
        clwl, clwr, clh, cdwl, cdwr, cdh, cdb = get_coefficients(
            aoaw + self.imain, aoaw + self.imain, aoah + self.itail
        )
        Lwl = 0.5 * self.rho * self.area[0] * clwl * v[0] ** 2
        Lwr = 0.5 * self.rho * self.area[1] * clwr * v[1] ** 2
        Lh = 0.5 * self.rho * self.area[2] * clh * v[2] ** 2
        Dwl = 0.5 * self.rho * self.area[0] * cdwl * v[0] ** 2
        Dwr = 0.5 * self.rho * self.area[1] * cdwr * v[1] ** 2
        Dh = 0.5 * self.rho * self.area[2] * cdh * v[2] ** 2
        Db = 0.5 * self.rho * self.area[3] * cdb * v[0] ** 2

        # Dwlz = 0.5 * self.rho * self.area[0] * cdwlz * self.velocity_linear[2] ** 2 * -np.sign(self.velocity_linear[2])
        # Dwrz = 0.5 * self.rho * self.area[1] * cdwrz * self.velocity_linear[2] ** 2 * -np.sign(self.velocity_linear[2])
        # Dhz = 0.5 * self.rho * self.area[2] * cdhz * self.velocity_linear[2] ** 2 * -np.sign(self.velocity_linear[2])

        wl = np.array([Dwlx, 0, -Lwl + Dwlz])
        wr = np.array([Dwrx, 0, -Lwr + Dwrz])
        h = np.array([Dhx, 0, -Lh + Dhz])

        wl = np.matmul(self.aero_to_body(pitch, beta), wl)
        wr = np.matmul(self.aero_to_body(pitch, beta), wr)
        h = np.matmul(self.aero_to_body(pitch, beta), h)
        b = np.matmul(self.aero_to_body(pitch, beta), b)

        return wl, wr, h, b

    def aero_to_body(self, aoa, b):
        """
        Args:
            aoa: angle of attack
            b: angle of sideslip
        Returns:
            T: the transformation matrix
        """
        a = aoa
        b = b
        T = np.array(
            [
                [np.cos(b) * np.cos(a), np.sin(b), np.cos(b) * np.sin(a)],
                [-np.sin(b) * np.cos(a), np.cos(b), -np.sin(b) * np.sin(a)],
                [-np.sin(a), 0, np.cos(a)],
            ]
        )
        T = np.linalg.inv(T)

        return T

    def moments(self, R, F):
        """
        Args:
            R: list of the distances for left wing, right wing, horizontal stabilizer, aileron left, aileron right, rudder and cg
            F: list of forces for the same
        Returns:
            M: sum of moment around the X axis, Y axis and Z axis
        """

        M = np.zeros((len(R), 3))
        for i in range(len(R)):
            M[i] = np.cross(R[i], F[i])
        M = np.array([sum(M.T[0]), sum(M.T[1]), sum(M.T[2])])
        Fnet = np.array([sum(F.T[0]), sum(F.T[1]), sum(F.T[2])])
        return M, Fnet

    def accelerations(self, F, M):
        """
        Args:
             F: force vector
             M: moment vector
        Returns:
            a: acceleration vector
            ang: angular acceleration vector
        """

        inv = np.linalg.inv(np.copy(self.I))
        temp = np.matmul(self.I, self.velocity_angular)
        rest = M - np.cross(self.velocity_angular, temp)
        dang_dt = np.matmul(inv, rest)
        temp2 = np.cross(self.velocity_angular, self.velocity_linear)
        a = np.array([sum(F.T[0]) / self.m, sum(F.T[1]) / self.m, sum(F.T[2]) / self.m]) - temp2
        return a, dang_dt

    def body_to_inertial(self):
        """
        Args:
        Returns:
            T: the transformation matrix
        """

        theta, gamma, psi = self.euler[1], self.euler[0], self.euler[2]

        T = np.array(
            [
                [np.cos(theta) * np.cos(psi), np.cos(theta) * np.sin(psi), -np.sin(theta)],
                [
                    np.sin(gamma) * np.sin(theta) * np.cos(psi) - np.cos(gamma) * np.sin(psi),
                    np.sin(gamma) * np.sin(theta) * np.sin(psi) + np.cos(gamma) * np.cos(psi),
                    np.sin(gamma) * np.sin(theta),
                ],
                [
                    np.cos(gamma) * np.sin(theta) * np.cos(gamma) + np.sin(gamma) * np.sin(psi),
                    np.cos(gamma) * np.sin(theta) * np.sin(psi) - np.sin(gamma) * np.cos(psi),
                    np.cos(gamma) * np.cos(theta),
                ],
            ]
        )

        T = np.linalg.inv(T)
        return T


class PID:
    def __init__(self, Kp=1.0, Ki=1.0, Kd=1.0):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd

        self.prev_error = 0.0
        self.integral = 0.0

    def __call__(self, error, dt):
        # Calculate the integral term
        self.integral += error * dt

        # Calculate the derivative term
        derivative = (error - self.prev_error) / dt

        # Update the previous error for the next iteration
        self.prev_error = error

        return self.Kp * error + self.Ki * self.integral + self.Kd * derivative


class Controller:
    def __init__(self, dt=0.01, Kp=0.0, Ki=0.0, Kd=0.0):
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
    def __init__(self, dt=0.01, Kp=0.0, Ki=0.0, Kd=0.0):
        self.PID_thrust_pitch = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, pitch):
        T = self.PID_thrust_pitch(pitch, self.dt)
        return T / 2, T / 2


class Controller3:
    def __init__(self, dt=0.01, Kp=0.0, Ki=0.0, Kd=0.0):
        self.PID_xthrust = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, xthrust):
        ithrust = self.PID_xthrust(xthrust, self.dt)
        return ithrust / 2, ithrust / 2


class Controller35:
    def __init__(self, dt=0.01, Kp=0.0, Ki=0.0, Kd=0.0):
        self.PID_xthrust = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, xthrust):
        ithrust = self.PID_xthrust(xthrust, self.dt)
        if abs(ithrust) > abs(system.max_thrust):
            ithrust = abs(system.max_thrust)
        return ithrust / 2, ithrust / 2


class Controller4:
    def __init__(self, dt=0.01, Kp=1.0, Ki=1.0, Kd=1.0):
        self.PID_theta = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_theta):
        ele = self.PID_theta(error_theta, self.dt)
        return ele


class Controller5:
    def __init__(self, dt=0.01, Kp=0.0, Ki=0.0, Kd=0.0):
        self.PID_thrust = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_mag_thrust, T_current, error_velocity_linear):
        T_change = self.PID_thrust(error_mag_thrust, self.dt)
        if abs(T_change) > abs(system.max_thrust):
            T_change = system.max_thrust

        T = np.array([0, 0, 0])
        T[2] = -(system.F[9][2] + system.F[0][2] + system.F[1][2] + system.F[2][2])
        controller_xthrust = Controller35(dt=dt, Kp=9.0, Ki=9.0, Kd=6.0)
        Tlx, Trx = controller_xthrust(error_velocity_linear[0])
        if abs(error_velocity_linear[0]) < 19:
            T[0] = Tlx + Trx
        else:
            T[0] = np.sqrt(T_change**2 - T[2] ** 2)
        return T


dt = 0.01
system = System()
controller_thrust = Controller(dt=dt, Kp=2000.0, Ki=0.0, Kd=300.0)
# controller_pitch = Controller2(dt=dt, Kp=0., Ki=0., Kd=0.)
controller_xthrust = Controller3(dt=dt, Kp=6000.0, Ki=2500.0, Kd=3000.0)
controller_elevator = Controller4(dt=dt, Kp=12000.0, Ki=600.0, Kd=4000.0)
controller_control_thrust = Controller5(dt=dt, Kp=1000.0, Ki=100.0, Kd=2500.0)

euler_ref = np.array([0.0, 0.0, 0.0])
velocity_linear_ref = np.array([400 / 3.6, 0.0, 0.0])
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

    Tl = np.array([Tlx, 0, Tlz])
    Tr = np.array([Trx, 0, Trz])
    T = Tl + Tr

    error_mag_thrust = abs(control_thrust_ref) - np.linalg.norm(T)

    T = controller_control_thrust(error_mag_thrust, T, error_velocity_linear)

    thrust_in_time.append(np.copy(T))

    exc_sig = T / 2, T / 2, Fele

    system(exc_sig, dt)

euler_in_time = np.array(euler_in_time)
velocity_linear_in_time = np.array(velocity_linear_in_time)
position_in_time = np.array(position_in_time)
thrust_in_time = np.array(thrust_in_time)

plt.plot(np.array(np.arange(0, 1e4, 1)) * dt, velocity_linear_in_time[:, 0])
plt.ylabel("Velocity in forward direction [m/s]")
plt.xlabel("Time [s]")
format_plot()
plt.show()
# plt.plot(velocity_linear_in_time[:, 1])
# plt.ylabel("velocity in the y direction")
# plt.show()
plt.plot(np.array(np.arange(0, 1e4, 1)) * dt, velocity_linear_in_time[:, 2])
plt.ylabel("velocity in the z direction")
plt.xlabel("Time [s]")
format_plot()
plt.show()

# plt.plot(position_in_time[:, 2], label="z-position")
# plt.plot(position_in_time[:, 0], label="x-position")
# plt.ylabel("position")
# plt.legend()
# plt.show()

plt.plot(
    np.array(np.arange(0, 1e4, 1)) * dt, thrust_in_time[:, 2], label="Thrust, downward positive"
)
plt.plot(
    np.array(np.arange(0, 1e4, 1)) * dt, thrust_in_time[:, 0], label="Thrust, forward positive"
)
plt.ylabel("Thrust [N]")
plt.xlabel("Time [s]")
plt.legend()
format_plot()
plt.show()

plt.plot(np.array(np.arange(0, 1e4, 1)) * dt, np.degrees(ithrust_in_time[:]))
plt.ylabel("Engine tilt angle [deg]")
plt.xlabel("Time [s]")
format_plot()
plt.show()

# plt.plot(euler_in_time[:, 2], color='b')
# plt.plot(euler_in_time[:, 1], color='g')
# plt.plot(euler_in_time[:, 0], color='r')
# plt.legend('ypr')
# plt.show()
