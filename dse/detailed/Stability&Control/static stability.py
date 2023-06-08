import numpy as np
import matplotlib.pyplot as plt
# import matplotlib as mpl
# mpl.use("TkAgg")


class Coefficients:

    def __init__(self):
        # for scissor plot
        self.mass = np.array([2350, 100, 250])  # cg, payload, pilots
        self.location = np.array([0.34, 0.94, -0.24])  # expressed in x/c
        # for tail sizing
        self.CL_alpha_h = 2 * np.pi  # cl alpha curve for tail
        self.CL_alpha_A = 2 * np.pi  # cl alpha curve for aircraft
        self.downwash_angle = 0  # downwash induced angle
        self.length_h = 11  # xh - xw (distance between tail and main wing)
        self.main_wing_chord = 4  # main wing chord
        self.Vh_V = 1  # ratio of tail speed to wing speed
        self.X_ac = 0.25  # location of aerodynamic center with respect to the chord
        self.SM = 0.05  # safety margin
        self.CL_h = -0.8  # -1 for full moving tail (what is it? lift of tail)
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
                         ((self.CL_alpha_h / self.CL_alpha_A) * (1 - self.downwash_angle) * (
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
                    (self.CL_alpha_h / self.CL_alpha_A) * (1 - self.downwash_angle) * (
                        self.length_h / self.main_wing_chord) * self.Vh_V ** 2)
        sh_s_control_value = (min(cgs) + self.Cmac / self.CL_A_h - self.X_ac) / (self.CL_h / self.CL_A_h) * (
                    self.length_h / self.main_wing_chord) * self.Vh_V ** 2
        sh_s_value = max(sh_s_control_value, sh_s_stability_value)
        cgrange = [min(cgs), max(cgs)]
        # print(sh_s_value, cgrange, "sh_s_value, cgrange")
        return sh_s_value, cgrange

    def neutral_stick_fixed(self):
        C_N_h_alpha = coefficients.CL_A_h
        C_N_alpha = coefficients.CL_alpha_A
        depsilon_dalpha = coefficients.downwash_angle
        V_h_V_ratio = coefficients.Vh_V
        Sh_S_ratio = coefficients.tail_area(cg)[0]
        l_h = coefficients.length_h
        c = coefficients.main_wing_chord
        x_w = coefficients.X_ac
        x_n_fix = (C_N_h_alpha / C_N_alpha) * (1 - depsilon_dalpha) * (V_h_V_ratio) ** 2 * (Sh_S_ratio * l_h / c) + x_w
        return x_n_fix

    def neutral_stick_free(self):
        C_N_h_alpha = coefficients.CL_A_h
        C_N_h_delta = 3.37  # Dummy number, TBD
        C_h_alpha = 0.00001  # Dummy number, TBD
        C_h_delta = -0.464  # Dummy number, TBD
        C_N_alpha = coefficients.CL_alpha_A
        depsilon_dalpha = coefficients.downwash_angle
        V_h_V_ratio = coefficients.Vh_V
        Sh_S_ratio = coefficients.tail_area(cg)[0]
        l_h = coefficients.length_h
        c = coefficients.main_wing_chord
        x_w = coefficients.X_ac
        C_N_h_alpha_free = C_N_h_alpha - C_N_h_delta * (C_h_alpha / C_h_delta)
        x_n_free = (C_N_h_alpha_free / C_N_alpha) * (1 - depsilon_dalpha) * (V_h_V_ratio) ** 2 * (Sh_S_ratio * l_h / c) + x_w
        return x_n_free


coefficients = Coefficients()
cg = coefficients.loading_diagram()
area_t, cg_range = coefficients.tail_area(cg)
neutral_point_fixed = coefficients.neutral_stick_fixed()
neutral_point_free = coefficients.neutral_stick_free()


class Equilibrium:

    def __init__(self):
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
        self.C_N_h = 1.0        # Normal coefficient of the tail
        self.C_M_ac_w = coefficients.Cmac   # Main wing moment coefficient
        self.C_M_ac_h = -0.2    # Tail moment coefficient
        # make sure that these are the same as the ones expressed in x/c under coefficients
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
        z_h = np.arange(-3, 0.1, 0.5)   # we pick z location
        wing = self.C_M_ac_w + self.C_N_w * self.X_W - self.C_T_w * self.Z_W
        thrust = (-Tz * self.X_T + Tx * self.Z_T) / (0.5 * self.rho * self.S * self.vel**2)
        tail = (self.C_M_ac_h * self.Chord_h - self.C_T_h * z_h) * ((self.Vh_V**2) * self.Sh_S / self.Chord_w)
        x_h = -(wing + thrust + tail) / (self.C_N_h * ((self.Vh_V**2) * self.Sh_S / self.Chord_w))
        # this is the minimum distance, when the tail is producing the maximum lift
        # if distance is larger, it just produces less lift
        return x_h[2]

    def rudder_sizing(self, x_rudder):
        v_cross_wind = 113 / 3.6
        s_fuselage_side = 15 * 3
        drag_coefficient_fuselage = 1.28
        safety_factor = 1.5
        side_slip_angle = np.arctan2(v_cross_wind, self.vel)
        cl_alpha = 2 * np.pi    # or just get the exact number at the side slip angle
        cl_due_to_sideslip = side_slip_angle * cl_alpha
        force_engine = 430 / 2
        distance_ratio = 17 / 10
        cl_deflection = 0.0504 * 180 / np.pi    # found in xflr5 will validate and verify the number
        max_deflection = np.radians(20)
        cl_due_to_deflection = cl_deflection * max_deflection
        D_fuselage = 0.5 * self.rho * s_fuselage_side * drag_coefficient_fuselage * v_cross_wind**2
        s_side_slip = (D_fuselage * safety_factor) / \
                      (0.5 * self.rho * self.vel**2 * (cl_due_to_sideslip + cl_due_to_deflection))
        s_engine = (force_engine * distance_ratio * safety_factor) / \
                   (0.5 * self.rho * self.vel**2 * cl_due_to_deflection)

        return s_side_slip, s_engine


equilibrium = Equilibrium()
thrust_x = equilibrium.sum_in_x()
thrust_z = equilibrium.sum_in_z()
x_location_tail = equilibrium.sum_moments(thrust_x, thrust_z)
s1, s2 = equilibrium.rudder_sizing(x_location_tail)
print(s1, s2)
