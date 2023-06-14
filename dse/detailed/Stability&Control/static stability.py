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
        self.CL_alpha_A = 5.271  # cl alpha curve for aircraft
        self.downwash_angle = 0  # downwash induced angle
        self.length_h = 13  # xh - xw (distance between tail and main wing)
        self.main_wing_chord = 3.11  # main wing chord
        self.Vh_V = 1  # ratio of tail speed to wing speed
        self.X_ac = 0.25  # location of aerodynamic center with respect to the chord
        self.SM = 0.05  # safety margin
        self.CL_h = -0.95  # -1 for full moving tail (what is it? lift of tail)
        self.CL_A_h = 1.33  # Cl of the aircraft with no tail?
        self.Cmac = -0.33  # main wing moment coefficient
        self.C_N_h_delta = 2.43     # Dummy number, TBD
        self.C_h_alpha = -0.118    # Dummy number, TBD
        self.C_h_delta = -0.279     # Dummy number, TBD
        self.C_h_delta_t = -0.228   # Dummy number, TBD
        self.C_m_0 = 0.08333325813508141    # Dummy number, TBD
        self.d_deltae_d_deltase = 2.18    # Dummy number, TBD
        self.W0 = 3000 * 3.71   # weight
        self.rho = 0.01         # density
        self.S = 133.5      # main wing area
        self.vel = 400 / 3.6    # velocity
        self.tail_chord = 1.897     # tail chord

        # # reference plane parameters
        # # for scissor plot
        # self.mass = np.array([1, 10000])  # cg, payload, pilots
        # self.location = np.array([0.409, 0.426])  # just min max cgs
        # self.CL_alpha_h = 5  # done
        # self.CL_alpha_A = 5.143  # done
        # self.downwash_angle = 0.25  # done
        # self.length_h = 4.03  # done
        # self.main_wing_chord = 1.47  # done
        # self.Vh_V = 0.9  # ratio of tail speed to wing speed
        # self.X_ac = 0.25  # done
        # self.SM = 0.1  # safety margin
        # self.CL_h = -0.35 * 2.83**(1/3)  # done
        # self.CL_A_h = 0.3  # done
        # self.Cmac = -0.04  # not sure
        # self.C_N_h_delta = 3.37  # Dummy number, TBD
        # self.C_h_alpha = 0.00001  # Dummy number, TBD
        # self.C_h_delta = -0.124  # Dummy number, TBD
        # self.C_h_delta_t = -0.124  # Dummy number, TBD
        # self.C_m_delta = -2  # Dummy number, TBD
        # self.C_m_0 = 0.1119  # done
        # self.d_deltae_d_deltase = 2  # Dummy number, TBD
        # self.W0 = 1043*9.81  # done
        # self.rho = 1.225 * 0.742 # done
        # self.S = 16.2  # done
        # self.vel = 226 / 3.6  # done
        # self.tail_chord = 0.707  # done

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
        plt.title("The loading diagram")
        plt.show()
        assert max(cgs) <= max(self.location), f"The aft cg is begind the most aft component"
        assert min(cgs) >= min(self.location), f"The forwards cg is in front of the foremost component"
        return cgs

    def tail_area(self, cgs):

        cgrange = np.arange(0, 1, 0.01)
        sh_s_stability = (cgrange - (self.X_ac - self.SM)) /\
                         ((self.CL_alpha_h / self.CL_alpha_A) * (1 - self.downwash_angle) * (
                            self.length_h / self.main_wing_chord) * self.Vh_V ** 2)
        sh_s_control = (cgrange + self.Cmac / self.CL_A_h - self.X_ac) / ((self.CL_h / self.CL_A_h) * (
                    self.length_h / self.main_wing_chord) * self.Vh_V ** 2)
        plt.plot(cgrange, sh_s_stability, color='tab:blue', label='Stability curve')
        plt.plot(cgrange, sh_s_control, color='tab:orange', label='Control curve')
        plt.axvline(x=max(cg), color='g', label='cg range')
        plt.axvline(x=min(cg), color='g')
        plt.ylim(0, 0.5)
        plt.legend()
        plt.show()

        sh_s_stability_value = (max(cgs) - (self.X_ac - self.SM)) / (
                    (self.CL_alpha_h / self.CL_alpha_A) * (1 - self.downwash_angle) * (
                        self.length_h / self.main_wing_chord) * self.Vh_V ** 2)
        sh_s_control_value = (min(cgs) + self.Cmac / self.CL_A_h - self.X_ac) / ((self.CL_h / self.CL_A_h) * (
                    self.length_h / self.main_wing_chord) * self.Vh_V ** 2)
        sh_s_value = max(sh_s_control_value, sh_s_stability_value)
        cgrange = [min(cgs), max(cgs)]
        print("sh_s_value:", sh_s_value, "cg range with respect to the chord", cgrange)
        assert sh_s_value > 0, f"The tail can not have negative volume, check the signs of the coefficients"
        assert sh_s_value < 1, f"The tail should be smaller than the main wing, check the given coefficients"
        return sh_s_value, cgrange

    def neutral_stick_fixed(self, areat):
        C_N_h_alpha = self.CL_alpha_h
        C_N_alpha = self.CL_alpha_A
        depsilon_dalpha = self.downwash_angle
        V_h_V_ratio = self.Vh_V
        Sh_S_ratio = areat
        l_h = self.length_h
        c = self.main_wing_chord
        x_w = self.X_ac
        x_n_fix = (C_N_h_alpha / C_N_alpha) * (1 - depsilon_dalpha) * V_h_V_ratio ** 2 * (Sh_S_ratio * l_h / c) + x_w
        print(f"x n fix is: {x_n_fix}")
        return x_n_fix

    def neutral_stick_free(self, areat):
        C_N_h_alpha = self.CL_alpha_h
        C_N_h_delta = self.C_N_h_delta  # Dummy number, TBD
        C_h_alpha = self.C_h_alpha  # Dummy number, TBD
        C_h_delta = self.C_h_delta  # Dummy number, TBD
        C_N_alpha = self.CL_alpha_A
        depsilon_dalpha = self.downwash_angle
        V_h_V_ratio = self.Vh_V
        Sh_S_ratio = areat
        l_h = self.length_h
        c = self.main_wing_chord
        x_w = self.X_ac
        C_N_h_alpha_free = C_N_h_alpha - C_N_h_delta * (C_h_alpha / C_h_delta)
        x_n_free = (C_N_h_alpha_free / C_N_alpha) * \
                   (1 - depsilon_dalpha) * (V_h_V_ratio ** 2) * (Sh_S_ratio * l_h / c) + x_w
        print(f"x n free is: {x_n_free}")
        return x_n_free

    def elevator_force(self, cg_chord, xn_free_chord, xn_fixed_chord, area_tail):
        C_m_delta = - self.C_N_h_delta * self.Vh_V**2 * area_tail * self.length_h / self.main_wing_chord
        single_deflection = (2 * self.W0 * self.C_h_delta * (cg_chord - xn_free_chord)) / \
                            (self.rho * self.S * (self.vel**2) * C_m_delta * self.C_h_delta)
        print(f"the elevator deflection is {np.degrees(single_deflection)}")
        V = np.arange(60, 150, 1)
        dFe_dV = self.d_deltae_d_deltase * area_tail * self.S * self.tail_chord * (self.Vh_V**2) * \
                 (-self.rho * self.vel * self.C_h_delta_t * single_deflection)
        print(f"the control force stability (dFe/dV) at Fe=0 is {dFe_dV} Ns/m")

        Fe = self.d_deltae_d_deltase * area_tail * self.S * self.tail_chord * self.Vh_V**2 * \
             (self.W0 * self.C_h_delta * (cg_chord - xn_free_chord) / (self.S * C_m_delta) - 0.5 * self.rho * V**2 * self.C_h_delta_t * single_deflection)

        plt.plot(V, Fe, color='tab:blue', label='Elevator force')
        plt.axvline(x=self.vel, color='tab:orange', label='cruise velocity')
        plt.title("elevator control force curve")
        plt.legend()
        plt.show()
        control_force = self.d_deltae_d_deltase * area_tail * self.S * self.tail_chord * self.Vh_V**2 * \
                        (self.W0 * self.C_h_delta * cg_chord * xn_free_chord / (self.S * C_m_delta) - 0.5 * self.rho * self.vel**2 * self.C_h_delta_t * single_deflection)

        elevator_deflection = (- 1 / C_m_delta) * (self.C_m_0 + self.W0 * (cg_chord - xn_fixed_chord) / (0.5 * self.rho * V**2 * self.S))
        plt.plot(V, elevator_deflection, color='tab:blue', label='Elevator trim curve')
        plt.axvline(x=self.vel, color='tab:orange', label='cruise velocity')
        plt.title("elevator trim curve")
        plt.legend()
        plt.show()

        return control_force


coefficients = Coefficients()
cg = coefficients.loading_diagram()
area_t, cg_range = coefficients.tail_area(cg)
neutral_point_fixed = coefficients.neutral_stick_fixed(area_t)
neutral_point_free = coefficients.neutral_stick_free(area_t)
elevator_force = coefficients.elevator_force(cg[0], neutral_point_free, neutral_point_fixed, area_t)

class Equilibrium:

    def __init__(self):
        self.C_T_w = 0.054      # tangential force coefficient of wing done
        self.C_T_h = 0.054      # tangential force coefficient of tail
        self.C_T_b = 0.1        # tangential force coefficient of body 0.029
        self.S = coefficients.S         # area of wing
        self.Sb = 40            # area of the body
        self.Sh_S = area_t     # area of the tail to wing ratio
        self.rho = coefficients.rho     # density
        self.vel = coefficients.vel     # velocity
        self.Vh_V = coefficients.Vh_V   # ratio of tail speed to wing speed
        self.W0 = coefficients.W0       # The weight of the aircraft
        self.pitch = 0          # The pitch angle of aircraft in radians done
        self.C_N_w = 1.33       # Normal coefficient of the wing (CL) at 0.6 degrees done
        self.C_N_h = 0.327        # Normal coefficient of the tail
        self.C_M_ac_w = coefficients.Cmac   # Main wing moment coefficient
        self.C_M_ac_h = -0.69   # Tail moment coefficient
        self.X_W = (max(cg) - coefficients.X_ac) * coefficients.main_wing_chord  # X location of main wing ac
        self.Z_W = -1           # Z location of the wing center of pressure (or aerodynamic center)
        self.X_T = 2            # X location of the thrust
        self.Z_T = -1           # Z location of the thrust
        self.Chord_w = coefficients.main_wing_chord     # Main wing chord
        self.Chord_h = coefficients.tail_chord        # Tail chord

        # # Verification and Validation: Cessna parameters
        # self.C_T_w = 0.02  # tangential force coefficient of wing done
        # self.C_T_h = 0.03  # tangential force coefficient of tail
        # self.C_T_b = 0.029  # tangential force coefficient of body
        # self.S = coefficients.S  # area of wing
        # self.Sb = 5.75  # area of the body
        # self.Sh_S = area_t + 0.03  # area of the tail to wing ratio
        # self.rho = coefficients.rho  # density
        # self.vel = coefficients.vel  # velocity
        # self.Vh_V = coefficients.Vh_V  # ratio of tail speed to wing speed
        # self.W0 = coefficients.W0  # The weight of the aircraft
        # self.pitch = 0  # The pitch angle of aircraft in radians done
        # self.C_N_w = 0.327  # aoa = 0.5 degrees for naca 2412
        # self.C_N_h = 0.156  # aoa = 1.2 degrees for naca 0012
        # self.C_M_ac_w = coefficients.Cmac  # Main wing moment coefficient
        # self.C_M_ac_h = -0.060  # Tail moment coefficient
        # # make sure that these are the same as the ones expressed in x/c under coefficients
        # self.X_W = (max(cg) - coefficients.X_ac) * coefficients.main_wing_chord  # X location of main wing ac
        # self.Z_W = -1  # Z location of the wing center of pressure (or aerodynamic center)
        # self.X_T = 2  # X location of the thrust
        # self.Z_T = 0  # Z location of the thrust
        # self.Chord_w = coefficients.main_wing_chord  # Main wing chord
        # self.Chord_h = coefficients.tail_chord  # Tail chord

    def sum_in_x(self):
        tx = (self.C_T_b * (self.Sb / self.S) + self.C_T_h * self.Sh_S * self.Vh_V**2 + self.C_T_w)\
             * (0.5 * self.rho * self.S * self.vel**2) + self.W0 * np.sin(self.pitch)
        print(f"tx is {tx}")
        return tx

    def sum_in_z(self):
        tz = (self.C_N_w + self.C_N_h * self.Sh_S * self.Vh_V**2)\
             * (0.5 * self.rho * self.S * self.vel**2) - self.W0 * np.cos(self.pitch)
        print(f"tz is {tz}")
        return tz

    def sum_moments(self, Tx, Tz):
        # we chose position of wing and engine and get position of the tail
        z_h = -1   # we pick z location
        print(f"main wing location {self.X_W}")
        Tz = 0
        wing = self.C_M_ac_w + self.C_N_w * self.X_W - self.C_T_w * self.Z_W
        thrust = (-Tz * self.X_T + Tx * self.Z_T) / (0.5 * self.rho * self.S * self.vel**2)
        tail = (self.C_M_ac_h * self.Chord_h - self.C_T_h * z_h) * ((self.Vh_V**2) * self.Sh_S / self.Chord_w)
        x_h = -(wing + thrust + tail) / (self.C_N_h * ((self.Vh_V**2) * self.Sh_S / self.Chord_w))

        momentwing = self.C_M_ac_w + self.C_N_w * self.X_W - self.C_T_w * self.Z_W
        momentthrust = (-Tz * self.X_T + Tx * self.Z_T) / (0.5 * self.rho * self.S * self.vel**2)
        momenttail = (self.C_M_ac_h * self.Chord_h - self.C_T_h * z_h + self.C_N_h * x_h) * ((self.Vh_V**2) * self.Sh_S / self.Chord_w)
        summoment = momenttail + momentthrust + momentwing
        print(f"sum of moments is {summoment}")
        return x_h, z_h

    def moments_change_alpha(self, x_h, z_h, alpha):
        C_N_h = self.C_N_h + coefficients.CL_alpha_h * alpha
        C_N_w = self.C_N_w + coefficients.CL_alpha_A * alpha
        C_T_w = self.C_T_w + 0.2962 * alpha
        C_T_h = self.C_T_h + 0.17425 * alpha

        Tz = (C_N_w + C_N_h * self.Sh_S * self.Vh_V ** 2) \
             * (0.5 * self.rho * self.S * self.vel ** 2) - self.W0 * np.cos(self.pitch)

        Tx = (self.C_T_b * (self.Sb / self.S) + C_T_h * self.Sh_S * self.Vh_V ** 2 + C_T_w) \
             * (0.5 * self.rho * self.S * self.vel ** 2) + self.W0 * np.sin(self.pitch)

        momentwing = self.C_M_ac_w + C_N_w * self.X_W - C_T_w * self.Z_W
        momentthrust = (-Tz * self.X_T + Tx * self.Z_T) / (0.5 * self.rho * self.S * self.vel ** 2)
        momenttail = (self.C_M_ac_h * self.Chord_h - C_T_h * z_h + C_N_h * x_h) * (
                    (self.Vh_V ** 2) * self.Sh_S / self.Chord_w)
        summoment = momenttail + momentthrust + momentwing
        # print(f"moment coefficient  : {summoment}")
        return summoment


    def rudder_sizing(self, x_rudder):
        # v_cross_wind = 55 / 3.6    # 7.7m for/s cesna
        # s_fuselage_side = 15 * 3
        # drag_coefficient_fuselage = 1.28
        # safety_factor = 1.5
        # side_slip_angle = np.arctan2(v_cross_wind, self.vel)
        # cl_alpha = 2 * np.pi
        # cl_due_to_sideslip = side_slip_angle * cl_alpha  # or just get the exact number at the sideslip angle
        # force_engine = 430 / 2
        # x_rudder = 10 * (x_rudder / x_rudder)   # remove later just so the number is not crazy
        # distance_ratio = (43.73/2) / x_rudder    # just 10 for now
        # cl_deflection = 0.0504 * 180 / np.pi    # found in xflr5 will validate and verify the number
        # max_deflection = np.radians(20)
        # cl_due_to_deflection = cl_deflection * max_deflection
        # d_fuselage = 0.5 * self.rho * s_fuselage_side * drag_coefficient_fuselage * v_cross_wind**2
        # s_side_slip = (d_fuselage * safety_factor) / \
        #               (0.5 * self.rho * self.vel**2 * (cl_due_to_deflection - cl_due_to_sideslip))
        # s_engine = (force_engine * distance_ratio * safety_factor) / \
        #            (0.5 * self.rho * self.vel**2 * cl_due_to_deflection)

        # Verification and Validation
        v_cross_wind = 7.7
        s_fuselage_side = 9.43
        drag_coefficient_fuselage = 1.28
        safety_factor = 1.
        side_slip_angle = np.arctan2(v_cross_wind, self.vel)
        cl_alpha = 5.7116
        cl_due_to_sideslip = side_slip_angle * cl_alpha
        cl_deflection = 0.0504 * 180 / np.pi
        max_deflection = np.radians(20)
        cl_due_to_deflection = cl_deflection * max_deflection
        d_fuselage = 0.5 * self.rho * s_fuselage_side * drag_coefficient_fuselage * v_cross_wind ** 2

        s_side_slip = (d_fuselage * safety_factor) / \
                      (0.5 * self.rho * self.vel ** 2 * (cl_due_to_deflection - cl_due_to_sideslip))

        s_engine = 0    # since only one engine
        return s_side_slip, s_engine


equilibrium = Equilibrium()
thrust_x = equilibrium.sum_in_x()
thrust_z = equilibrium.sum_in_z()
x_location_tail, z_location_tail = equilibrium.sum_moments(thrust_x, thrust_z)
print("minimum x location of tail", x_location_tail)
s1, s2 = equilibrium.rudder_sizing(x_location_tail)
print(f"minimum tail size for sideslip {s1}, and engine off {s2}")

cmintime = []
for aoa in np.arange(np.radians(0), np.radians(5), np.radians(0.1)):
    a = equilibrium.moments_change_alpha(x_location_tail, z_location_tail, aoa)
    cmintime.append(a)

plt.plot(np.arange(0, 5, 0.1), cmintime)
plt.title("Cm vs alpha")
plt.show()
def landing_distance():
    m = 3000
    g = 3.71
    thrust_up = 5000
    rho = 0.01
    surface_area_w = 135
    cl_max = 2.4
    surface_area_b = 20
    cd_body = 0.1   # get from pedro
    mu_break = 0.5  # performance brakes coefficient
    v_stall = np.sqrt((m * g - thrust_up) / (0.5 * rho * surface_area_w * cl_max))
    v_approach = v_stall * 1.1
    v = v_approach
    distance = 0
    dt = 0.01

    while v > 0:
        drag = (m * g - thrust_up) * mu_break + 0.5 * rho * cd_body * surface_area_b * v**2
        a = drag/m
        v = v - a * dt
        distance += v * dt

        if thrust_up > 0:
            thrust_up = thrust_up - 5

    print("landing distance", distance)  # first order estimate ignoring the cl of the wing


landing_distance()

def wing_stuff(area, aspect, taper):
    b = np.sqrt(area * aspect)
    c = area/b
    temp = c * 2 / (taper+1)
    croot = temp * taper
    ctip = temp
    print(f"b = {b}, croot = {croot}, ctip = {ctip}")

print("Main wing:")
wing_stuff(coefficients.S, 15, 2)
print("Horizontal tail:")
wing_stuff(coefficients.S*(area_t + equilibrium.Sh_S), 7.5, 2.5)
print("Vertical tail:")
wing_stuff(11.33, 2, 2.5)
print(coefficients.S)
print(coefficients.S*(area_t + equilibrium.Sh_S))
