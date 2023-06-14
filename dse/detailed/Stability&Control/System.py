import numpy as np
import Aircraft


class System:
    def __init__(self):
        self.geometry = Aircraft.define_geometry()
        self.area, self.imain, self.itail, self.m, self.I, self.W0 = Aircraft.define_areas()
        # self.rho = 0.01   # for Mars
        self.rho = 1.225 * 0.742    # for Cessna

        # Test for verification and validation Cesssna 172
        self.euler_prev = np.array([0, np.radians(10), 0])  # the euler angles: roll(at X), pitch(at Y), yaw(at Z)
        self.velocity_linear_prev = np.array([226/3.6, 0, 0.])
        self.velocity_angular_prev = np.array([0, 0, 0.])

        # Test for stability our plane
        # self.euler_prev = np.array([0, 0, 0])  # the euler angles: roll(at X), pitch(at Y), yaw(at Z)
        # self.velocity_linear_prev = np.array([400/3.6, 0, 0.])
        # self.velocity_angular_prev = np.array([0, np.radians(15), 0.])

        # Test for controllers
        # self.euler_prev = np.array([np.radians(5), np.radians(-5), np.radians(5)])  # the euler angles: roll(at X), pitch(at Y), yaw(at Z)
        # self.velocity_linear_prev = np.array([105, 0, 0.])
        # self.velocity_angular_prev = np.array([0.1, 0, 0.])

        # current state
        self.euler = np.copy(self.euler_prev)
        self.velocity_linear = np.copy(self.velocity_linear_prev)
        self.velocity_angular = np.copy(self.velocity_angular_prev)

    def get_state(self):
        return self.euler, self.velocity_linear, self.velocity_angular

    def __call__(self, exc_sig, dt):
        Fal, Far, Fe, Fr, Tlx, Trx = exc_sig
        aoawings = np.arctan2(self.velocity_linear[2] - self.velocity_angular[1], self.velocity_linear[0])
        aoatail = np.arctan2(self.velocity_linear[2] + self.velocity_angular[1]*10, self.velocity_linear[0])
        velocitywings = np.array([self.velocity_linear[0] - self.velocity_angular[1], self.velocity_linear[0] - self.velocity_angular[1], self.velocity_linear[0] + self.velocity_angular_prev[1] * 10])
        wl, wr, h = self.aero_forces(aoawings, aoatail, velocitywings, self.euler[2], self.euler[1])
        W = np.array([-self.W0 * np.sin(self.euler[1]), 0, self.W0 * np.cos(self.euler[1])])
        F = np.array([wl, wr, h, np.array(Fal), np.array(Far), np.array([0, 0, Fe]), np.array(Fr), np.array([Tlx, 0, 0]),
                      np.array([Trx, 0, 0]), W])
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
        clwl, clwr, clh, cdwl, cdwr, cdh = Aircraft.get_coefficients(aoaw + self.imain, aoaw + self.imain, aoah + self.itail)
        Lwl = 0.5 * self.rho * self.area[0] * clwl * v[0] ** 2
        Lwr = 0.5 * self.rho * self.area[1] * clwr * v[1] ** 2
        Lh = 0.5 * self.rho * self.area[2] * clh * v[2] ** 2
        Dwl = 0.5 * self.rho * self.area[0] * cdwl * v[0] ** 2
        Dwr = 0.5 * self.rho * self.area[1] * cdwr * v[1] ** 2
        Dh = 0.5 * self.rho * self.area[2] * cdh * v[2] ** 2

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
        # a = aoa
        # b = b
        # T = np.array([[np.cos(b) * np.cos(a), np.sin(b), np.cos(b) * np.sin(a)],
        #               [-np.sin(b) * np.cos(a), np.cos(b), -np.sin(b) * np.sin(a)],
        #               [-np.sin(a), 0, np.cos(a)]])
        # T = np.linalg.inv(T)
        T = np.array([[np.cos(a) * np.cos(b), -np.cos(a) * np.sin(b), -np.sin(a)],
                     [np.sin(b), np.cos(b), 0],
                     [np.sin(a) * np.cos(b), -np.sin(a) * np.sin(b), np.cos(a)]])

        return T

    def moments(self, R, F):
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
        Fnet = np.array([sum(F.T[0]), sum(F.T[1]), sum(F.T[2])])

        #Rwingleft, Rwingright, Rstabilizer, Raileronleft, Raileronright, Relevator, Rrudder, Rthrustleft, Rthrustright, Cg])
        return M, Fnet

    def accelerations(self, F, M):
        '''

        Args:
             F: force vector
             M: moment vector
        Returns:
            a: acceleration vecotor
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


class SecondOrderSystem:
    def __init__(self):
        self.m = 1.
        self.c = 0.1
        self.k = 0.1
        self.A = np.array([[0., 1.],
                           [-self.k/self.m, -self.c/self.m]])
        self.B = np.array([[0.],
                           [1/self.m]])
        self.x = np.array([1., -0.5])

    @property
    def F_c(self):
        # return np.array([np.cos(self.x[0])])
        return np.array([0.78])

    def step(self, F_u, dt):
        dx_dt = np.matmul(self.A, self.x) + np.matmul(self.B, F_u + self.F_c)
        self.x = self.x + dx_dt*dt
        return self.x