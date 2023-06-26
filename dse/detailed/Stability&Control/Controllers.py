# import numpy as np


class PID:
    def __init__(self, Kp=1.0, Ki=1.0, Kd=1.0):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd

        self.prev_error = 0.0
        self.integral = 0.0

    def __call__(self, error, dt):
        # Calculate the integral term
        self.integral += error * dt

        # Calculate the derivative term
        derivative = (error - self.prev_error) / dt

        # Update the previous error for the next iteration
        self.prev_error = error

        return self.Kp * error + self.Ki * self.integral + self.Kd * derivative


class Controller:
    def __init__(self, dt=0.01, Kp=1.0, Ki=1.0, Kd=1.0):
        self.PID_theta = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_theta):
        ele = self.PID_theta(error_theta, self.dt)
        if ele > 1277:
            ele = 1277
        elif ele < -1277:
            ele = -1277
        return ele


class Controller2:
    def __init__(self, dt=0.01, Kp=1.0, Ki=1.0, Kd=1.0):
        self.PID_thrust = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_velocity):
        T = self.PID_thrust(error_velocity, self.dt)
        # for our plane
        Tlx = T / 2
        Trx = T / 2
        if Tlx > 900:
            Tlx = 900
            Trx = 900
        elif Tlx < 0:
            Tlx = 0
            Trx = 0
        return Tlx, Trx

        # # for Cessna 172 verification and validation
        # return T


class Controller3:
    def __init__(self, dt=0.01, Kp=1.0, Ki=1.0, Kd=1.0):
        self.PID_rudder = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_yaw):
        rud = self.PID_rudder(error_yaw, self.dt)
        if rud > 731:
            rud = 731
        elif rud < -731:
            rud = -731
        return rud


class Controller4:
    def __init__(self, dt=0.01, Kp=1.0, Ki=1.0, Kd=1.0):
        self.PID_aileron = PID(Kp=Kp, Ki=Ki, Kd=Kd)
        self.dt = dt

    def __call__(self, error_roll):
        ail = self.PID_aileron(error_roll, self.dt)
        if ail > 1066:
            ail = 1066
        elif ail < -1066:
            ail = -1066
        return ail / 2, -ail / 2
