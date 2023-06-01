from System import System
from Controllers import Controller, Controller2, Controller3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("TkAgg")

# state_now
# reference_now
# send error to controller
# get controller's control signal
# send control signal to system
# update new state

dt = 0.01
system = System()
controller_theta = Controller(dt=dt, Kp=25000., Ki=6000., Kd=500.)
controller_thrust = Controller2(dt=dt, Kp=1000., Ki=5., Kd=0.)
controller_phi = Controller3(dt=dt, Kp=-1000., Ki=5., Kd=-5000.)

euler_ref = np.array([0, 0, 0.])
velocity_linear_ref = np.array([400/3.6, 0, 0.])
euler_in_time = []
velocity_linear_in_time = []
#Tlx, Trx = 200, 200
# velprev = 110

for i in range(int(1e4)):
    euler, velocity_linear, velocity_angular = system.get_state()

    euler_in_time.append(np.copy(euler))
    velocity_linear_in_time.append(np.copy(velocity_linear))

    error_euler = euler_ref - euler
    error_velocity_linear = velocity_linear_ref - velocity_linear

    Tlx, Trx = controller_thrust(error_velocity_linear[0]) #think about if we control just x velocity or everything
    if Tlx > 430:
        Tlx = 430
        Trx = 430
    rud = controller_phi(error_euler[2])
    print(rud)
    elev = controller_theta(error_euler[1])
    exc_sig = np.array([0, 0, 0.]), np.array([0, 0, 0.]), elev, np.array([0, rud, 0.]), Tlx, Trx

    #print(controller_theta(error_euler[1]), error_euler[1])

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
plt.plot(euler_in_time[:, 2], color='b')
plt.plot(euler_in_time[:, 1], color='g')
plt.plot(euler_in_time[:, 0], color='r')
plt.legend('ypr')
plt.show()

