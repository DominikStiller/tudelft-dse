import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

df = pd.read_csv(filepath_or_buffer="lift_distribution.csv", 
                 sep=",",
                 skipinitialspace=True,
                 skiprows=20,
                 header=0,
                 nrows=38,
                 )

b = 15.
rho = 0.01
V_inf = 111
del_V_inf = 5.
N = 1000
y_interpolated = np.linspace(-b, b, N)
dy_interpolated = y_interpolated[1] - y_interpolated[0]
V_interpolated = np.full(N, V_inf)
V_interpolated[(y_interpolated < -b + 10.4) | (y_interpolated > b - 10.4)] = V_inf + del_V_inf

plt.plot(y_interpolated, V_interpolated)
plt.show()

y_span = df["y-span"].to_numpy()
Chord = df["Chord"].to_numpy()
Cl = df["Cl"].to_numpy()
CL = pd.read_csv(filepath_or_buffer="lift_distribution.csv", 
                 sep=",",
                 skipinitialspace=True,
                 skiprows=9,
                 header=None,
                 nrows=1,
                 ).iloc[0, 1]

y_boundaries = [-b]
for i, y in enumerate(y_span):
    y_boundaries.append(y_boundaries[i] + 2*(y-y_boundaries[i]))

y_boundaries = np.array(y_boundaries)
dy = y_boundaries[1:] - y_boundaries[:-1]

# Cl_interpolated = np.interp(y_interpolated, y_span, Cl)
cs_Cl = CubicSpline(y_span, Cl)
Cl_interpolated = cs_Cl(y_interpolated)
cs_Chord = CubicSpline(y_span, Chord)
Chord_interpolated = cs_Chord(y_interpolated)

# plt.plot(y_span, Cl)
# plt.plot(y_interpolated, Cl_interpolated)
# plt.show()

C_L_check_interpolated = sum(Cl_interpolated*dy_interpolated/(2*b))
corrective_constant = CL / C_L_check_interpolated

L_interpolated = sum(dy_interpolated*Chord_interpolated*0.5*rho*V_interpolated**2*Cl_interpolated) * corrective_constant

print(L_interpolated, corrective_constant)

L_needed = 3000 * 3.71
print(L_needed / L_interpolated)


