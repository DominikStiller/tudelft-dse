import sys
sys.path.append("../../dse")

from plotting import format_plot, save_plot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### thin airfoil theory graph start ###

df = pd.read_csv("T1_Re100.000_M0.10_N9.0.csv", delimiter=",", header=0)
alpha = df.iloc[:, 0].to_numpy()
cl = df.iloc[:, 1].to_numpy()

thin_airfoil_theory = 2 * np.pi**2 * alpha / 180

error_absolute = thin_airfoil_theory - cl
error_relative = error_absolute / thin_airfoil_theory

print("absolute maximum error", max(error_absolute))
print("relative maximum error", max(error_relative))

fig, ax = plt.subplots(1, 1, figsize=(8, 3), sharey="all")


plt.plot(alpha, thin_airfoil_theory, label="thin airfoil theory - analytical solution", linewidth=3)
plt.plot(alpha, cl, label="NACA0006 - XFLR5 simulation", linewidth=3)
plt.xlabel("angle of attack [deg]")
plt.ylabel("coefficient of lift [-]")
plt.legend()

format_plot()
save_plot(".", "thin_airfoil_theory")

plt.show()

### thin airfoil theory graph end ###


### XFLR5 S1123 graph start ###



### XFLR5 S1123 graph end ###

