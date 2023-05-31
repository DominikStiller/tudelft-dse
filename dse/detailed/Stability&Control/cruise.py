from System import System
from Controllers import Controller, thrust_Controller
import numpy as np
import matplotlib.pyplot as plt


# state_now
# reference_now
# send error to controller
# get controller's control signal
# send control signal to system
# update new state

dt = 0.01
system = System()
controller_theta = Controller(dt=dt, Kp=1., Ki=0., Kd=0.)

euler_ref = np.array([0, 0, 0.])
velocity_linear_ref = np.array([111, 0, 0.])
euler_in_time = []
velocity_linear_in_time = []

for i in range(1000):
    euler, velocity_linear, velocity_angular = system.get_state()

    euler_in_time.append(np.copy(euler))
    velocity_linear_in_time.append(np.copy(velocity_linear))

    error_euler = euler_ref - euler
    error_velocity_linear = velocity_linear_ref - velocity_linear

    Tlx, Trx = thrust_Controller(-1, -1, -1, -1, -1, -1.)
    exc_sig = np.array([0, 0, 0.]), np.array([0, 0, 0.]), controller_theta(error_euler[1]), np.array([0, 0, 0.]), Tlx, Trx

    system(exc_sig, dt)

euler_in_time = np.array(euler_in_time)
velocity_linear_in_time = np.array(velocity_linear_in_time)

plt.plot(velocity_linear_in_time[:, 0])
plt.show()
