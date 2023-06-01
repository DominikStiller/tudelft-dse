import numpy as np


class PID:
    def __init__(self, Kp=1., Ki=1., Kd=1.):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd

        self.prev_error = 0.
        self.integral = 0.

    def __call__(self, error, dt):

        # Calculate the integral term
        self.integral += error * dt

        # Calculate the derivative term
        derivative = (error - self.prev_error) / dt

        # Update the previous error for the next iteration
        self.prev_error = error

        return self.Kp * error + self.Ki * self.integral + self.Kd * derivative


class Controller:
    def __init__(self, dt=0.01, Kp=1., Ki=1., Kd=1.):
        self.PID_theta = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_theta):
        return self.PID_theta(error_theta, self.dt)

class Controller2:
    def __init__(self, dt=0.01, Kp=1., Ki=1., Kd=1.):
        self.PID_thrust = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_velocity):
        T = self.PID_thrust(error_velocity, self.dt)
        return T/2, T/2




