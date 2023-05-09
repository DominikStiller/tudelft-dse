import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit


def eliminate_nan(df, x_key, y_key):
    df = df[df[x_key].notna()]
    df = df[df[y_key].notna()]
    return df


# Read the excel and store data in a dataframe
similar_vehicles = pd.read_excel('Tiltrotor_examples.xlsx', sheet_name='Data', index_col=0)
similar_vehicles = similar_vehicles.loc[:, similar_vehicles.columns != 'source']

# Convert units
similar_vehicles['shaft power'] = similar_vehicles['shaft power'] / 1.341

# Divide into manned and unmanned examples
manned_vehicles = similar_vehicles.loc[similar_vehicles['crew_size'] > 0]
unmanned_vehicles = similar_vehicles.loc[similar_vehicles['crew_size'] == 0]

### Power - Propeller diameter ###
# df1 = eliminate_nan(manned_vehicles, 'propeller diameter', 'shaft power')
# df2 = eliminate_nan(unmanned_vehicles, 'propeller diameter', 'shaft power')
#
# # plt.scatter(df1['propeller diameter'], df1['shaft power'], label='manned')
# # plt.scatter(df2['propeller diameter'], df2['shaft power'], label='unmanned')
# # plt.xlabel('Propeller diameter [m]')
# # plt.ylabel('Shaft power [kW]')
# # plt.legend()
# # plt.show()
#
# # Obtain a relationship a*e^(b*d)+c
# diameter_power =pd.concat([df1, df2])
# popt, pcov = curve_fit(lambda d, a, b, c: a * np.exp(b * d) + c, diameter_power['propeller diameter'], diameter_power['shaft power'])
#
# print(f'std = {np.sqrt(np.diag(pcov))}')
#
# a, b, c = popt
# diameter_fitted = np.linspace(np.min(diameter_power['propeller diameter']),
#                               np.max(diameter_power['propeller diameter']), 100)
# power_fitted = a * np.exp(b * diameter_fitted) + c
# plt.scatter(diameter_power['propeller diameter'], diameter_power['shaft power'], color='r')
# plt.plot(diameter_fitted, power_fitted, label=f'P = {round(a, 2)}exp({round(b, 2)}*D)+({round(c, 2)})')
# plt.xlabel('Propeller diameter [m]')
# plt.ylabel('Shaft power [kW]')
# plt.legend()
# plt.show()


### useful mass - MTOM ###
df1 = eliminate_nan(manned_vehicles, 'useful mass', 'MTOM')
df2 = eliminate_nan(unmanned_vehicles, 'useful mass', 'MTOM')
combined_dataframe = pd.concat([df1, df2])
#
# # plt.scatter(df1['useful mass'], df1['MTOM'], label='manned')
# # plt.scatter(df2['useful mass'], df2['MTOM'], label='unmanned')
# # plt.xlabel('Useful mass [kg]')
# # plt.ylabel('MTOM [kg]')
# # plt.legend()
# # plt.show()
#
# # Find regression
# regression = linregress(combined_dataframe['useful mass'], combined_dataframe['MTOM'])
#
# useful_mass_fitted = np.linspace(np.min(combined_dataframe['useful mass']),
#                                  np.max(combined_dataframe['useful mass']), 100)
# MTOM_fitted = regression.slope * useful_mass_fitted + regression.intercept
#
# plt.scatter(combined_dataframe['useful mass'], combined_dataframe['MTOM'], color='r')
# plt.plot(useful_mass_fitted, MTOM_fitted, label='MTOM = ' + f'{round(regression.slope, 2)}' + r'm$_u$ + ' + f'{round(regression.intercept, 2)}')
# plt.xlabel('Useful mass [kg]')
# plt.ylabel('MTOM [kg]')
# plt.legend()
# plt.show()


### MTOM - Power ###

# df1 = eliminate_nan(manned_vehicles, 'shaft power', 'MTOM')
# df2 = eliminate_nan(unmanned_vehicles, 'shaft power', 'MTOM')
# combined_dataframe = pd.concat([df1, df2])

# plt.scatter(df1['MTOM'], df1['shaft power'], label='manned')
# plt.scatter(df2['MTOM'], df2['shaft power'], label='unmanned')
# plt.xlabel('MTOM [kg]')
# plt.ylabel('Shaft power [kW]')
# plt.legend()
# plt.show()

# regression = linregress(combined_dataframe['MTOM'], combined_dataframe['shaft power'])
#
# MTOM_fitted = np.linspace(np.min(combined_dataframe['MTOM']),
#                           np.max(combined_dataframe['MTOM']), 100)
# power_fitted = regression.slope * MTOM_fitted + regression.intercept
#
# plt.scatter(combined_dataframe['MTOM'], combined_dataframe['shaft power'], color='r')
# plt.plot(MTOM_fitted, power_fitted, label='P = ' + f'{round(regression.slope, 2)}' + 'MTOM + ' + f'({round(regression.intercept, 2)})')
# plt.xlabel('Useful mass [kg]')
# plt.ylabel('MTOM [kg]')
# plt.legend()
# plt.show()