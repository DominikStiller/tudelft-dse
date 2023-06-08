import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("TkAgg")


class Coefficients:

    def __init__(self):
        # for scissor plot
        self.mass = np.array([2350, 100, 250])
        self.location = np.array([0.34, 0.94, -0.24])  # expressed in x/c
        # for tail sizing
        self.CL_alph_ah = 1  # cl alpha curve for tail
        self.CL_alpha_A = 1  # cl alpha curve for aircraft
        self.downwash_angle = 0  # downwash induced angle
        self.length_h = 11  # xh - xw (distance between tail and main wing)
        self.main_wing_chord = 4  # main wing chord
        self.Vh_V = 1  # ratio of tail speed to wing speed
        self.X_ac = 0.2  # location of aerodynamic center with respect to the chord
        self.SM = 0.05  # safety margin
        self.CL_h = -1  # -1 for full moving tail (what is it? lift of tail)
        self.CL_A_h = 1.8  # Cl of the aircraft with no tail?
        self.Cmac = -0.2  # main wing moment coefficient

    def loading_diagram(self):
        cgs = np.zeros(len(self.mass))
        tot_mass = np.zeros(len(self.mass))
        for i in range(len(self.mass)):
            if tot_mass[i - 1] != 0:
                cgs[i] = (cgs[i - 1] * tot_mass[i - 1] + self.location[i] * self.mass[i]) / (
                            self.mass[i] + tot_mass[i - 1])
            else:
                cgs[i] = self.location[i]

            tot_mass[i] = tot_mass[i - 1] + self.mass[i]

        plt.plot(cgs, tot_mass)
        plt.xlabel('x/c')
        plt.ylabel('total mass [kg]')
        # plt.show()
        return cgs

    def tail_area(self, cgs):

        cgrange = np.arange(0, 1, 0.01)
        sh_s_stability = (cgrange - (self.X_ac - self.SM)) /\
                         ((self.CL_alph_ah / self.CL_alpha_A) * (1 - self.downwash_angle) * (
                            self.length_h / self.main_wing_chord) * self.Vh_V ** 2)
        sh_s_control = (cgrange + self.Cmac / self.CL_A_h - self.X_ac) / (self.CL_h / self.CL_A_h) * (
                    self.length_h / self.main_wing_chord) * self.Vh_V ** 2

        plt.plot(cgrange, sh_s_stability, color='tab:blue', label='Stability curve')
        plt.plot(cgrange, sh_s_control, color='tab:orange', label='Control curve')
        plt.axvline(x=max(cg), color='g', label='cg range')
        plt.axvline(x=min(cg), color='g')
        plt.ylim(0, 0.5)
        plt.legend()
        plt.plot()
        # plt.show()

        sh_s_stability_value = (max(cgs) - (self.X_ac - self.SM)) / (
                    (self.CL_alph_ah / self.CL_alpha_A) * (1 - self.downwash_angle) * (
                        self.length_h / self.main_wing_chord) * self.Vh_V ** 2)
        sh_s_control_value = (min(cgs) + self.Cmac / self.CL_A_h - self.X_ac) / (self.CL_h / self.CL_A_h) * (
                    self.length_h / self.main_wing_chord) * self.Vh_V ** 2
        sh_s_value = max(sh_s_control_value, sh_s_stability_value)
        cgrange = [min(cgs), max(cgs)]
        print(sh_s_value, cgrange, "sh_s_value, cgrange")
        return sh_s_value, cgrange


coefficients = Coefficients()
cg = coefficients.loading_diagram()
area_t, cg_range = coefficients.tail_area(cg)


class Equilibrium:

    def __init__(self):
        # self.Tx = 800           # total thrust in x direction
        self.C_T_w = 0.05       # tangential force coefficient of wing
        self.C_T_h = 0.05       # tangential force coefficient of tail
        self.C_T_b = 0.01       # tangential force coefficient of body
        self.S = 112            # area of wing
        self.Sb = 40            # area of the body
        self.Sh_S = area_t      # area of the tail to wing ratio
        self.rho = 0.01         # density
        self.vel = 400 / 3.6    # velocity
        self.Vh_V = coefficients.Vh_V   # ratio of tail speed to wing speed
        self.W0 = 2700 * 3.71   # The weight of the aircraft
        self.pitch = 0          # The pitch angle of aircraft in radians
        self.C_N_w = 1.3553        # Normal coefficient of the wing
        self.C_N_h = 1.2        # Normal coefficient of the tail
        self.C_M_ac_w = coefficients.Cmac   # Main wing moment coefficient
        self.C_M_ac_h = -0.2    # Tail moment coefficient
        self.X_W = 1            # X location of the wing center of pressure (or aerodynamic center) im not sure
        self.Z_W = -1           # Z location of the wing center of pressure (or aerodynamic center)
        self.X_T = 2            # X location of the thrust
        self.Z_T = -1           # Z location of the thrust
        self.Chord_w = coefficients.main_wing_chord     # Main wing chord
        self.Chord_h = 2        # Tail chord

    def sum_in_x(self):
        tx = (self.C_T_b * (self.Sb / self.S) + self.C_T_h * self.Sh_S * self.Vh_V**2 + self.C_T_w)\
             * (0.5 * self.rho * self.S * self.vel**2) + self.W0 * np.sin(self.pitch)
        return tx

    def sum_in_z(self):
        tz = (self.C_N_w + self.C_N_h * self.Sh_S * self.Vh_V**2)\
             * (0.5 * self.rho * self.S * self.vel**2) - self.W0 * np.cos(self.pitch)
        return tz

    def sum_moments(self, Tx, Tz):
        # we chose position of wing and engine and get position of the tail
        z_h = np.arange(0, 3.1, 0.5)
        wing = self.C_M_ac_w + self.C_N_w * self.X_W + self.C_T_w * self.Z_W
        thrust = (Tz * self.X_T - Tx * self.Z_T) / (0.5 * self.rho * self.S * self.vel**2)
        print(Tz, self.X_T, Tx, self.Z_T, self.rho, self.S, self.vel**2)
        tail = (self.C_M_ac_h * self.Chord_h + self.C_T_h * z_h) * ((self.Vh_V**2) * self.Sh_S / self.Chord_w)
        x_h = (wing + thrust + tail) / (self.C_N_h * ((self.Vh_V**2) * self.Sh_S / self.Chord_w))
        return x_h


equilibrium = Equilibrium()
thrust_x = equilibrium.sum_in_x()
thrust_z = equilibrium.sum_in_z()
M = equilibrium.sum_moments(thrust_x, thrust_z)
print(M)
