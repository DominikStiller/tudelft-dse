from System import System
from Controllers import Controller, Controller2, Controller3, Controller4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("TkAgg")

dt = 0.01
system = System()
controller_pitch = Controller(dt=dt, Kp=25000., Ki=6000., Kd=5000.)
controller_thrust = Controller2(dt=dt, Kp=1000., Ki=5., Kd=0.)
controller_yaw = Controller3(dt=dt, Kp=-15000., Ki=-50., Kd=-15000.)
controller_roll = Controller4(dt=dt, Kp=-4000., Ki=-1., Kd=-7500.)

euler_ref = np.array([0, 0, 0.])
velocity_linear_ref = np.array([400/3.6, 0, 0.])
euler_in_time = []
velocity_linear_in_time = []

# # for Cessna verification and validation
# velocity_linear_ref = np.array([226/3.6, 0, 0.])

for i in range(int(5e4)):
    euler, velocity_linear, velocity_angular = system.get_state()
    euler_in_time.append(np.copy(euler))
    velocity_linear_in_time.append(np.copy(velocity_linear))

    error_euler = euler_ref - euler
    error_velocity_linear = velocity_linear_ref - velocity_linear

    elev = controller_pitch(error_euler[1])
    # for Mars plane
    Tlx, Trx = controller_thrust(error_velocity_linear[0])

    # # for Cessna verification and validation
    # Tlx = controller_thrust(error_velocity_linear[0])
    # Trx = 0

    rud = controller_yaw(error_euler[2])
    Fal, Far = controller_roll(error_euler[0])

    # exc_sig = [0, 0, Fal], [0, 0, Far], elev, [0, rud, 0.], Tlx, Trx
    exc_sig = [0, 0, 0], [0, 0, 0], 0, [0, 0, 0], Tlx, Trx
    system(exc_sig, dt)

euler_in_time = np.array(euler_in_time)
velocity_linear_in_time = np.array(velocity_linear_in_time)

plt.plot(velocity_linear_in_time[:, 0])
plt.ylabel("velocity in the x direction")
plt.show()
plt.plot(velocity_linear_in_time[:, 2])
plt.ylabel("velocity in the z direction")
plt.show()
plt.plot(velocity_linear_in_time[:, 1])
plt.ylabel("velocity in the y direction")
plt.show()
plt.plot(np.degrees(euler_in_time[:, 2]), color='b')
plt.plot(np.degrees(euler_in_time[:, 1]), color='g')
plt.plot(np.degrees(euler_in_time[:, 0]), color='r')
plt.legend('ypr')
plt.show()
