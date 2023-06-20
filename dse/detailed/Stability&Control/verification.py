import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from System import SecondOrderSystem
from Controllers import PID

mpl.use("MacOSX")

system = SecondOrderSystem()
Ku = 10
Tu = 2
controller = PID(Kp=1, Ki=0.1, Kd=1)
dt = 0.01
steps = int(1e4)
x_in_time = [system.x]

x_ref = 4.0

for i in range(steps):
    x_in_time.append(system.step(controller(error=x_ref - system.x[0], dt=dt), dt))

x_in_time = np.array(x_in_time)

plt.plot(x_in_time[:, 0], label="position")
plt.plot(x_in_time[:, 1], label="velocity")
plt.plot(np.full(shape=steps, fill_value=system.F_c / system.k), label="equilibrium wo controller")
plt.plot(np.full(shape=steps, fill_value=x_ref), label="reference state")
plt.legend()
plt.show()
