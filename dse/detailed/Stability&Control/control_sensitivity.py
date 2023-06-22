import numpy as np
import pandas as pd
from SALib.analyze import sobol
from SALib.sample import saltelli
from tqdm import tqdm

from static_stability import Coefficients


def calculate_parameters(inputs):
    coeff = Coefficients(inputs)
    cgs = coeff.loading_diagram()
    shs, cgr = coeff.tail_area(cgs)
    xn_fix = coeff.neutral_stick_fixed(shs)
    xn_free = coeff.neutral_stick_free(shs)
    control_force = coeff.elevator_force(cgs[0], xn_free, xn_fix, shs)
    return shs, xn_fix, xn_free, control_force


if __name__ == "__main__":
    coefficients = [
        "mass_low",
        "location_low",
        "Cl_alpha_h",
        "Cl_alpha_a",
        "downwash_angle",
        "length_h",
        "main_wing_chord",
        "Vh_V",
        "X_ac",
        "SM",
        "CL_h",
        "Cl_Ah",
        "C_mac",
        "C_N_h_delta",
        "C_h_alpha",
        "C_h_delta",
        "C_h_delta_t",
        "C_m_delta",
        "Cm0",
        "d_deltae_d_deltase",
        "W0",
        "airDensity",
        "S",
        "Vel",
        "tailChord",
        "location_up",
        "mass_up",
    ]
    bounds = [
        [1, 5000],  # mass lower bound
        [0.4, 0.45],  # location lower bound
        [3, 7],  # CL_alpha_h
        [3, 7],  # CL_alpha_A
        [0, 0.99],  # downwash_angle
        [5, 15],  # length_h
        [1, 5],  # main_wing_chord
        [0.1, 0.99],  # Vh_V
        [0.1, 0.5],  # X_ac
        [0.05, 0.25],  # SM
        [-0.8, -0.1],  # CL_h
        [0.5, 2],  # CL_A_h
        [-2, -0.01],  # Cmac
        [1, 5],  # C_N_h_delta
        [1e-6, 1e-4],  # C_h_alpha
        [-0.5, -0.1],  # C_h_delta
        [-0.5, -0.1],  # C_h_delta_t
        [-5, -0.5],  # C_m_delta
        [0.1, 1.5],  # C_m_0
        [0.5, 5],  # d_deltae_d_deltase
        [800 * 9.81, 1200 * 9.81],  # W0
        [0.6, 1.225],  # rho
        [60, 200],  # S
        [75, 150],  # vel
        [1, 4],  # tail_chord
        [0.45, 0.5],  # Location upper bound
        [5000, 10000],  # Mass upper bound
    ]

    Outputs = ["Tail area", "xn_fix", "xn_free", "control_force"]

    problem = {"num_vars": len(bounds), "names": coefficients, "bounds": bounds}

    param_values = saltelli.sample(problem, 2**6)
    # Sh_S
    Y = np.zeros((param_values.shape[0], 4))
    for i in tqdm(range(len(param_values))):
        Y[i] = calculate_parameters(param_values[i])

    Si = [sobol.analyze(problem, y) for y in Y.T]

    S = np.zeros((len(Si), np.shape(Si[0]["ST"])[0]))
    S_err = np.zeros((len(Si), np.shape(Si[0]["ST"])[0]))
    for i in range(len(Si)):
        S[i] = Si[i]["ST"]
        S_err[i] = Si[i]["ST_conf"]

    df1 = pd.DataFrame(np.round(S.T, 2), index=coefficients, columns=Outputs)
    df1.to_excel("Coefficients_sensitivity.xlsx")
    print(df1.to_latex())
