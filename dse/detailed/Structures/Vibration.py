import numpy as np
import matplotlib.pyplot as plt
import vibration_toolbox as vtb
from scipy.integrate import odeint
import control
from material_properties import *

class TailVibes():
    def __init__(self):
        # Material properties
        self.E = materials['CFRP'].E
        self.rho = materials['CFRP'].rho
        # Beam geometry
        self.r = 0.1  # [m] Radius of the tail beam
        self.l = 10  # [m] length of the beam
        self.t_skin = 0.005  # [m] thickness of the beam skin
        self.m_beam = self.rho * np.pi * ((self.r+self.t_skin) ** 2 - (self.r-self.t_skin)**2) * self.l      # [kg] Mass of the tail beam
        self.I = np.pi * (2 * self.r)**3 * self.t_skin/8
        self.m_tail = 50



    def simsetting(self):
        # Simulation Parameters
        t_start = 0  # [s]
        t_end = 5  # [s]
        self.dt = 1e-6  # [s] Time step
        self.t = np.arange(t_start, t_end + self.dt, self.dt)

    def tailplanegeom(self, VH):
        Sv = 20  # The surface area of the vertical tail
        Sh = 20  # The surface area of the horizontal tail
        Cd = 1.28
        if VH == 'V':
            S = Sv
        elif VH == 'H':
            S = Sh
        else:
            print('The vibrations are not correctly designed')
            S = 1
        return S, Cd
    def sysparam(self,VH):
        S, Cd = self.tailplanegeom(VH)
        rho0 = 0.01
        V = 400/3.6
        q = 0.5 * rho0 * V ** 2
        self.ceq = (Cd * S * 0.5 * rho0) / self.m_tail
        self.keq = ((2 * np.pi * 3 / (2 * self.l) * q * S + (3 * self.E * self.I) / (self.l ** 3))) / self.m_tail  # Stiffness of the spring

    def userinput(self, VH, ah):
        S = self.tailplanegeom(VH)[0]
        q = 0.5 * 0.01 * (400/3.6) ** 2
        self.i = 3
        if self.i == 0:
            self.F_u = (-2 * np.pi * q * S * ah) / self.m_tail * np.ones(len(self.t))
        elif self.i== 1:
            self.F_u = (-2 * np.pi * q * S * ah) / self.m_tail * np.hstack((np.arange(0, 1 + self.dt, self.dt), np.ones(len(self.t)-len(np.arange(0, 1 + self.dt, self.dt)))))
        elif self.i == 2:
            self.F_u = np.zeros(len(self.t))
            self.F_u[0] = 1/self.dt
        elif self.i ==3:
            w = np.pi
            i_end = int(2/self.dt)
            self.F_u = np.sin(w * self.t[:i_end])
            self.F_u = np.hstack((self.F_u, np.zeros(len(self.t)-i_end)))
        print(f"The applied aerodynamic load is {max(abs(self.F_u)) * self.m_tail} [N]")

    def syssim(self):
        X_init = np.array([0., 0.])
        X = []
        X.append(X_init)

        for i in range(len(self.t) - 1):
            x_i = X[-1][0]
            x_dot_i = X[-1][1]
            X.append(X[-1] + self.dt * np.array([x_dot_i, self.F_u[i] - self.ceq * np.sign(x_dot_i) * x_dot_i ** 2 - self.keq * x_i]))
        self.x = np.array(X)[:, 0]
        self.v = np.array(X)[:, 1]

    def results(self):
        period = []
        for i in range(int(len(self.t)/2), len(self.x)-1):
            if self.x[i-1] < self.x[i] > self.x[i+1]:
                period.append(i)
                print('k')
        P = self.t[period[2]] - self.t[period[1]]
        self.avg = (max(self.x[int(len(self.t)/2):]) + min(self.x[int(len(self.t)/2):])) / 2
        deflection = max(self.x[int(len(self.t)/2):] - self.avg)
        print(f"The period of the response is {P} [s]")
        print(f"The natural frequency is {1/P} [Hz]")
        print(f"The extra deflection is {deflection} [m]")
        print(f"The new steady state is {self.avg} [m]")
        self.pi = 3 * self.E * self.I / self.l**3 * deflection
        print(f"The induced load due to the vibration is {self.pi} [N]")
    def plot(self):
        nth = int((1 / self.dt) / 1e4)
        plt.plot(self.t[::nth], self.x[::nth], label="Displacement")
        plt.plot(self.t[::nth], self.v[::nth], label="Velocity")
        # plt.plot(self.t[::nth], self.F_u[::nth]/self.keq, label="Force")
        # plt.axhline(self.avg)
        plt.legend()
        plt.show()



T = TailVibes()


S, Cd = T.tailplanegeom(VH='V')
T.simsetting()
T.sysparam(VH='V')
T.userinput(VH='V', ah=0.1)
T.syssim()
T.results()
T.plot()
print(S, Cd)


#
# print('k')
#
#
#
#
#
# # I = 0.5 * m_beam * ((r + t_skin)**2 - (r - t_skin)**2)            # [kgm2] Mass moment of inertia of the beam
#
# # Tail geometry
# Sh = 20  # [m2] The surface area of the tail
# m_tail = 50 - m_beam
# Cd = 1.28  # [-] The drag
#
# # Free stream properties
# rho0 = 0.01  # [kgm3] Airdensity
# V = 400/3.6  # [m/s] Cruise speed
# q = 0.5 * rho0 * V ** 2
#
#
#
#
#
#
# # Simulation Parameters
# t_start = 0  # [s]
# t_end = 10  # [s]
# dt = 1e-6  # [s] Time step
# t = np.arange(t_start, t_end + dt, dt)
#
# # Parameters defining the system
# m = m_tail + m_beam  # Mass
# ceq = (Cd * Sh * 0.5 * rho0)/m
# keq = ((2 * np.pi * 3/(2 * l) * q * Sh + (3 * E * I)/(l**3)))/m  # Stiffness of the spring
# # User input
# ah = 0.2 ##### User input
# F_u = (-2 * np.pi * q * Sh * ah)/m * np.ones(len(t))
#
#
#
#
# X_init = np.array([0., 0.])
# X = []
# X.append(X_init)
#
# for i in range(len(t)-1):
#     x_i = X[-1][0]
#     x_dot_i = X[-1][1]
#     X.append(X[-1] + dt * np.array([x_dot_i, F_u[i] - ceq * np.sign(x_dot_i) * x_dot_i**2 - keq * x_i]))
#
# X = np.array(X)
# nth = int((1/dt) / 1e4)
# plt.plot(t[::nth], X[:,0][::nth], label="Displacement")
# plt.plot(t[::nth], X[:,1][::nth], label="Velocity")
# plt.legend()
# plt.show()
