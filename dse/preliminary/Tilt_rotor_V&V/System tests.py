from dse.preliminary.Tilt_rotor.tiltrotormain import design
from SALib.sample import saltelli
from SALib.analyze import sobol
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import numpy as np


def test_function_sensitivity():
    outputs = ['Total mass', 'Rotor mass', 'Battery mass', 'Solar panel mass', 'Body mass', 'Wing mass', 'Tail mass',
               'Rotor radius', 'Wingspan', 'Chord', 'Cruise thrust', 'ROC_cruise', 'ROC_rotor']
    functions = ['Rotor sizing', 'Wing area sizing', 'Drag computation', 'Battery and panel sizing', 'Class2Weight']
    original_data = np.array(design(iterate=True))

    data = np.zeros((5, 13))
    for i in range(5):
        data[i] = 100 * (design(iterate=True, functionSensitivity=(i, 1.1)) - original_data) / original_data

    df = pd.DataFrame(np.round(data.T, 2), index=outputs, columns=functions)
    df.to_excel('F_sensitivity/Function_sensitivity.xlsx')
    print(df.to_latex())

    for i in tqdm(range(len(data.T))):
        plt.bar(functions, data.T[i])
        plt.title(outputs[i])
        plt.xticks(rotation=15, ha='right')
        plt.ylabel('% impact on final result')
        plt.tight_layout()
        plt.savefig(f'F_sensitivity/Function_sensitivity_{outputs[i]}.pdf')
        plt.show()


if __name__ == '__main__':
    test_function_sensitivity()

    bounds = [[0, 2], [0, 0.5], [0.01, 0.02], [100, 120], [200, 240], [400, 500], [1200, 1500], [0.05, 0.1],
              [200, 400], [300, 400], [0.8, 1]]

    outputs = ['Total mass', 'Rotor mass', 'Battery mass', 'Solar panel mass', 'Body mass', 'Wing mass', 'Tail mass',
               'Rotor radius', 'Wingspan', 'Chord', 'Cruise thrust', 'ROC_cruise', 'ROC_rotor']

    inputs = ['cl', 'cd', 'rho', 'V_cr', 'a', 'E_density', 'P_density', 'fill factor', 'TO_time', 'payload', 'oswald']

    problem = {
        'num_vars': 11,
        'names': ['cl', 'cd', 'airDensity', 'cruiseSpeed', 'soundSpeed', 'takeoffBatteryEnergyDensity',
                  'takeoffBatteryPowerDensity', 'fillFactor', 'takeOffTime', 'payloadMass', 'oswald'],
        'bounds': bounds
    }

    param_values = saltelli.sample(problem, 2**6)

    Y = np.zeros((param_values.shape[0], len(outputs)))
    print('Started computing sensitivities')
    for i in tqdm(range(len(param_values))):
        Y[i] = design(iterate=True, inputSensitivity=param_values[i])

    Si = [sobol.analyze(problem, y) for y in Y.T]
    print('Sensitivities computed')

    S = np.zeros((len(Si), np.shape(Si[0]['ST'])[0]))
    for i in range(len(Si)):
        S[i] = Si[i]['ST']

    df = pd.DataFrame(np.round(S.T, 2), index=inputs, columns=outputs)
    df.to_excel('I_sensitivity/Input_output_correlation.xlsx')
    print(df.to_latex())

    df = pd.DataFrame(np.round(S, 2), index=outputs, columns=inputs)
    df.to_excel('I_sensitivity/Input_output_correlation.xlsx')
    print(df.to_latex())

    for i in range(len(Si)):
        total_Si, first_Si, second_Si = Si[i].to_df()
        plt.bar(inputs, total_Si.ST.array)
        plt.xticks(rotation=15, ha='right')
        plt.title(outputs[i])
        plt.tight_layout()
        plt.savefig(f'I_sensitivity/Input_sensitivity_{outputs[i]}.pdf')
        plt.show()

