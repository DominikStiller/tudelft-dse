#import numpy as np


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


def thrust_Controller(v, aoa, a, Tl, Tr, Vcruise):
    return 160., 160.  # so that one day this is a PID

    # '''
    #
    # Args:
    #     v: velocity vector
    #     aoa: angle of attack vector including pitch at index 3
    #     a: acceleration vector
    #     Tl: thrust of the left endine
    #     Tr: thrust of the right engine
    # Returns:
    #     Tl: thrust of the left engine
    #     Tr: thrust of the right engine
    # '''
    # if a[0]>0.05:
    #     #accelerating too fast
    #     Tl = Tl - 50
    #     Tr = Tr - 50
    # elif a[0]<-0.05:
    #     #deaccelerating to fast
    #     Tr = Tr + 50
    #     Tl = Tl + 50
    # elif v[0] > Vcruise * np.cos(np.radians(aoa[1])) and a[0] > 0:
    #     #accelerating away from cruise
    #     Tr = Tr - 20
    #     Tl = Tl - 20
    # elif v[0] < Vcruise * np.cos(np.radians(aoa[1])) and a[0] < 0:
    #     #deaccelerating away from cruise
    #     Tr = Tr + 20
    #     Tl = Tl + 20
    # elif v[0] < Vcruise * np.cos(np.radians(aoa[1])) and a[0] > 0:
    #     #accelerating towards cruise
    #     Tr = Tr - 10 * abs(v[0] - Vcruise * np.cos(np.radians(aoa[1])))/abs(Vcruise * np.cos(np.radians(aoa[1])))
    #     Tl = Tl - 10 * abs(v[0] - Vcruise * np.cos(np.radians(aoa[1])))/abs(Vcruise * np.cos(np.radians(aoa[1])))
    # elif v[0] > Vcruise * np.cos(np.radians(aoa[1])) and a[0] < 0:
    #     #deacelerating towards cruise
    #     Tr = Tr + 10 * abs(v[0] - Vcruise * np.cos(np.radians(aoa[1]))) / abs(Vcruise * np.cos(np.radians(aoa[1])))
    #     Tl = Tl + 10 * abs(v[0] - Vcruise * np.cos(np.radians(aoa[1]))) / abs(Vcruise * np.cos(np.radians(aoa[1])))
    # else:
    #     Tl = Tl
    #     Tr = Tr
    #
    # # return Tl, Tr
