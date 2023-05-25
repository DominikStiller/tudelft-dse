import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from SALib.analyze import sobol
from SALib.sample import saltelli
from tqdm import tqdm

from dse.preliminary.Tilt_rotor.tiltrotormain import design


def test_function_sensitivity(Print=False):
    outputs = [
        "Range",
        "Total mass",
        "Rotor mass",
        "Battery mass",
        "Body mass",
        "Wing mass",
        "Tail mass",
        "Rotor radius",
        "Wingspan",
        "Chord",
        "Cruise thrust",
        "Max power",
    ]
    functions = [
        "Rotor sizing",
        "Wing area sizing",
        "Drag computation",
        "Battery and panel sizing",
        "Class2Weight",
    ]
    original_data = np.array(design(iterate=True))

    data = np.zeros((len(functions), len(outputs)))
    for i in range(5):
        data[i] = (
            100
            * (design(iterate=True, functionSensitivity=(i, 1.1)) - original_data)
            / original_data
        )

    df = pd.DataFrame(np.round(data.T, 2), index=outputs, columns=functions)
    print(df.to_latex())

    if Print:
        df.to_excel("F_sensitivity/Function_sensitivity.xlsx")
        for i in range(len(data)):
            plt.bar(outputs, data[i])
            plt.title(functions[i])
            plt.xticks(rotation=20, ha="right")
            plt.ylabel("% impact on final result")
            plt.tight_layout()
            plt.savefig(f"F_sensitivity/Function_sensitivity_{functions[i]}.pdf")
            plt.show()


if __name__ == "__main__":
    # test_function_sensitivity()

    outputs = [
        "Range",
        "Total mass",
        "Rotor mass",
        "Battery mass",
        "Body mass",
        "Wing mass",
        "Tail mass",
        "Rotor radius",
        "Wingspan",
        "Chord",
        "Cruise thrust",
        "Max power",
    ]

    # Define inputs to be analyzed
    givenValues = ["payload mass", "design range", "cruise speed", "maximum mass"]
    bounds1 = [[200, 350], [8e5, 1e6], [110, 150], [2700, 3500]]

    calculatedValues = ["gravity Mars", "air density", "speed of sound", "cruise viscosity"]
    bounds2 = [[3.6, 3.71], [0.01, 0.02], [215, 250], [4e-4, 6e-4]]

    assumedValues = [
        "cl",
        "cd",
        "t/c",
        "battery power density",
        "battery energy density",
        "blade density",
        "fill factor",
        "ultimate load",
        "oswald",
        "take off time",
    ]
    bounds3 = [
        [1.6, 2],
        [0.01, 0.1],
        [0.01, 0.02],
        [1200, 1500],
        [430, 600],
        [1300, 1525],
        [0.07, 0.085],
        [1, 1.75],
        [0.9, 1],
        [300, 600],
    ]

    # Analyze the sensitivity of the system to the requirements
    problem1 = {
        "num_vars": 4,
        "names": ["payloadMass", "designRange", "cruiseSpeed", "maxMass"],
        "bounds": bounds1,
    }

    param_values1 = saltelli.sample(problem1, 2**6)

    Y = np.zeros((param_values1.shape[0], len(outputs)))
    for i in tqdm(range(len(param_values1))):
        Y[i] = design(iterate=True, inputSensitivity=(problem1["names"], param_values1[i]))

    Si1 = [sobol.analyze(problem1, y) for y in Y.T]

    S = np.zeros((len(Si1), np.shape(Si1[0]["ST"])[0]))
    S_err = np.zeros((len(Si1), np.shape(Si1[0]["ST"])[0]))
    for i in range(len(Si1)):
        S[i] = Si1[i]["ST"]
        S_err[i] = Si1[i]["ST_conf"]

    df1 = pd.DataFrame(np.round(S.T, 2), index=givenValues, columns=outputs)
    # df1.to_excel('I_sensitivity/Given_values_sensitivity.xlsx')
    print(df1.to_latex())

    # Print results
    for i in range(len(S.T)):
        plt.bar(outputs, S.T[i], yerr=S_err.T[i], capsize=5)
        plt.xticks(rotation=15, ha="right")
        plt.title(givenValues[i])
        plt.tight_layout()
        plt.savefig(f"I_sensitivity/Customer_sensitivity_{givenValues[i]}.pdf")
        plt.show()

    # Analyze the sensitivity of the system to calculated parameters
    problem2 = {
        "num_vars": 4,
        "names": ["gravityMars", "airDensity", "soundSpeed", "visc_cr"],
        "bounds": bounds2,
    }

    param_values2 = saltelli.sample(problem2, 2**6)

    Y = np.zeros((param_values2.shape[0], len(outputs)))
    for i in tqdm(range(len(param_values2))):
        Y[i] = design(iterate=True, inputSensitivity=(problem2["names"], param_values2[i]))

    Si2 = [sobol.analyze(problem2, y) for y in Y.T]

    S = np.zeros((len(Si2), np.shape(Si2[0]["ST"])[0]))
    S_err = np.zeros((len(Si2), np.shape(Si2[0]["ST"])[0]))
    for i in range(len(Si2)):
        S[i] = Si2[i]["ST"]
        S_err[i] = Si2[i]["ST_conf"]

    df2 = pd.DataFrame(np.round(S.T, 2), index=calculatedValues, columns=outputs)
    df2.to_excel("I_sensitivity/Calculated_values_sensitivity.xlsx")
    print(df2.to_latex())

    # Print results
    for i in range(len(S.T)):
        plt.bar(outputs, S.T[i], yerr=S_err.T[i], capsize=5)
        plt.xticks(rotation=15, ha="right")
        plt.title(calculatedValues[i])
        plt.tight_layout()
        plt.savefig(f"I_sensitivity/Calculated_values_sensitivity_{calculatedValues[i]}.pdf")
        plt.show()

    # Analyze the sensitivity of the system to the assumed values
    problem3 = {
        "num_vars": 10,
        "names": [
            "cl",
            "cd",
            "t/c",
            "takeoffBatteryPowerDensity",
            "takeoffBatteryEnergyDensity",
            "bladeDensity",
            "fillFactor",
            "ultimateLoad",
            "takeOffTime",
            "oswald",
        ],
        "bounds": bounds3,
    }

    param_values3 = saltelli.sample(problem3, 2**6)

    Y = np.zeros((param_values3.shape[0], len(outputs)))
    for i in tqdm(range(len(param_values3))):
        Y[i] = design(iterate=True, inputSensitivity=(problem3["names"], param_values3[i]))

    Si3 = [sobol.analyze(problem3, y) for y in Y.T]

    S = np.zeros((len(Si3), np.shape(Si3[0]["ST"])[0]))
    S_err = np.zeros((len(Si3), np.shape(Si3[0]["ST"])[0]))
    for i in range(len(Si3)):
        S[i] = Si3[i]["ST"]
        S_err[i] = Si3[i]["ST_conf"]

    df3 = pd.DataFrame(np.round(S.T, 2), index=assumedValues, columns=outputs)
    df3.to_excel("I_sensitivity/Assumed_values_correlation.xlsx")
    print(df3.to_latex())

    # Print results
    for i in range(len(S.T)):
        plt.bar(outputs, S.T[i], yerr=S_err.T[i], capsize=5)
        plt.xticks(rotation=15, ha="right")
        plt.title(assumedValues[i])
        plt.tight_layout()
        if assumedValues[i] == "t/c":
            assumedValues[i] = "tc"
        plt.savefig(f"I_sensitivity/Assumed_values_sensitivity_{assumedValues[i]}.pdf")
        plt.show()
