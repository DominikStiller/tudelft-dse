import matplotlib.pyplot as plt
import numpy as np
import scipy



aoa = np.radians(5)
Fx = []
Dx = []
Fy = []
Dy = []
Fz = []
Dz = []
I = [20000, 20000, 20000] #moment of inertia around x, y, z axis
rho = 0.01
clalpha = 2 * np.pi
clalphah = 2 * np.pi
clcd = 20
clcdh = 20
S = [100, 10]
m = 3000
gmars = 3.71
Tmax = 12000
Vclimb = 2

#simulation
duration = 20 #duration of simulation in seconds
dt = 0.1 #time increment in seconds
Tangvel = np.radians(3*dt) #3 deg/s ang vel of thrust

def add_force(F,pos_vector):
    """

    :param F: Perturbing force vector in array form
    :param pos_vector: Position of force w.r.t. center of gravity in array form
    :return: Force and position vector
    """

def thrust_force(mode, Faero, W, ang):
    if mode == 0:
        T = np.array([0, 0, -Tmax])
        return T
    elif mode == 1:
        T = -1* (W + Faero)
        return T
    elif mode == 2:
        T = -1 * W - np.array([0, 0, Faero[2]])
        newang = ang + Tangvel
        Told = T[0]
        T[0] = T[2] / np.tan(newang)
        if abs(Tmax) > np.linalg.norm(T):
            return T, newang
        else:
            T[0] = Told
            return T, newang

    elif mode == 3:
        Tcurrent = Faero[0] + m * ang
        T = [Tcurrent, 0, 0]
        return T




def aerodynamic_force(V, alpha, beta, alphah, S):
    CLw = clalpha * alpha
    CDw = CLw / clcd

    CLh = clalphah * alphah
    CDh = CLh / clcdh

    Cdwz = 0.2
    Cdhz = 0.2

    q = 0.5 * rho

    if V[2] < 0:
        Fwz = - CLw * q * S[0] * V[0] ** 2 + (Cdwz * q * S[0] * V[2] ** 2)
        Fhz = - CLh * q * S[1] * V[0] ** 2 + (Cdhz * q * S[1] * V[2] ** 2)
    else:
        Fwz = - CLw * q * S[0] * V[0] ** 2 - (Cdwz * q * S[0] * V[2] ** 2)
        Fhz = - CLh * q * S[1] * V[0] ** 2 + (Cdhz * q * S[1] * V[2] ** 2)

    Fwx = - CDw * q * S[0] * V[0] ** 2
    Fhx = - CDh * q * S[1] * V[0] ** 2

    Fw = np.array([Fwx, 0, Fwz]).T
    Fh = np.array([Fhx, 0, Fhz]).T


    Fw = np.matmul(aero_to_body(alpha, beta), Fw).T
    Fh = np.matmul(aero_to_body(alpha, beta), Fh).T

    return Fw, Fh

def moments_Y(Fy, Dy):
    x = 2




def angular_acceleration(M, I):
    """
    Calculate the angular acceleration of aircraft.

    Args:
        I: Moments of inertia around X, Y, Z axis [kgm^2]
        M: Moment around X, Y, Z axis [rad]

    Returns:
        angular acceleration around X, Y, Z axis [rad/s^2]
    """

    alphaX = M[0]/I[0]
    alphaY = M[1]/I[1]
    alphaZ = M[2]/I[2]
    angacc = [alphaX, alphaY, alphaZ]

    return angacc


def accelerations(F, m):
    """
    Calculate the angular acceleration of aircraft.

    Args:
        m: mass of aircraft [kg]
        F: Forces X, Y, Z axis [N]

    Returns:
        acceleration of Z, Y, Z axis [m/s^2]
    """
    sumX = np.array([])
    sumY = np.array([])
    sumZ = np.array([])
    for i in range(len(F)):
        sumX.append(F[i][0])
        sumY.append(F[i][1])
        sumZ.append(F[i][2])

    ax = sumX/m
    ay = sumY/m
    az = sumZ/m
    acc = [ax, ay, az]
    return acc

def aero_to_body(a, b):
    T = np.array([[np.cos(b) * np.cos(a), np.sin(b), np.cos(b) * np.sin(a)],
         [-np.sin(b) * np.cos(a), np.cos(b), -np.sin(b) * np.sin(a)],
         [-np.sin(a), 0, np.cos(a)]])
    T = np.linalg.inv(T)

    return T

def body_to_inertial(theta, gamma, psi):
    T = np.array([[np.cos(theta) * np.cos(psi), np.cos(theta) * np.sin(psi), -np.sin(theta)],
                  [np.sin(gamma) * np.sin(theta) * np.cos(psi) - np.cos(gamma) * np.sin(psi), np.sin(gamma) * np.sin(theta) * np.sin(psi) + np.cos(gamma) * np.cos(psi), np.sin(gamma) * np.sin(theta)],
                  [np.cos(gamma) * np.sin(theta) * np.cos(gamma) + np.sin(gamma) * np.sin(psi), np.cos(gamma) * np.sin(theta) * np.sin(psi) - np.sin(gamma) * np.cos(psi), np.cos(gamma) * np.cos(theta)]])

    T = np.linalg.inv(T)

    return T

def run_simulation(duration, dt):
    pitch = 0
    alpha = 0    # will be a function of pitch and the velocity in body
    alphah = alpha # will be a function of alpha
    beta = 0
    W0 = m * gmars
    t = np.zeros(int(duration / dt))
    V = np.zeros((int(duration / dt), 3))
    X = np.zeros((int(duration / dt), 3))
    A = np.zeros((int(duration / dt), 3))
    Vi = np.array([0, 0, 0])
    Xi = np.array([0, 0, 0])
    Ai = np.array([0, 0, 0])
    ti = 0
    V[0] = Vi
    X[0] = Xi
    A[0] = Ai
    t[0] = ti
    i = 0
    # initial acceleration phase

    while V[i][2] > -Vclimb:
        i = i + 1
        mode = 0
        W = np.array([-W0 * np.sin(pitch), 0, W0 * np.cos(pitch)])
        Fw, Fh = aerodynamic_force(V[i - 1], alpha, beta, alphah, S)
        T = thrust_force(mode, (Fw + Fh), W, 0)
        Fnet = W + Fw + Fh + T
        A[i] = Fnet / m
        V[i] = V[i - 1] + A[i] * dt
        X[i] = X[i - 1] + V[i] * dt
        t[i] = np.round(t[i - 1] + dt, 2)
        if t[i] >= duration - dt:
            return A, V, X, t


    # climb until 500m phase
    while X[i][2] > -500:
        i = i + 1
        mode = 1
        W = np.array([-W0 * np.sin(pitch), 0, W0 * np.cos(pitch)])
        Fw, Fh = aerodynamic_force(V[i - 1], alpha, beta, alphah, S)
        T = thrust_force(mode, Fw+Fh, W, 0)
        Fnet = W + Fw + Fh + T
        A[i] = Fnet / m
        V[i] = V[i - 1] + A[i] * dt
        X[i] = X[i - 1] + V[i] * dt
        t[i] = np.round(t[i - 1] + dt, 2)
        if t[i] >= duration - dt:
            return A, V, X, t

    V[i][2] = 0
    Tang = np.radians(-90)
    while Tang<0:
        i = i + 1
        mode = 2

        W = np.array([-W0 * np.sin(pitch), 0, W0 * np.cos(pitch)])
        Fw, Fh = aerodynamic_force(V[i - 1], alpha, beta, alphah, S)
        T, Tang = thrust_force(mode, Fw + Fh, W, Tang)
        Fnet = W + Fw + Fh + T
        A[i] = Fnet / m
        V[i] = V[i - 1] + A[i] * dt
        X[i] = X[i - 1] + V[i] * dt
        t[i] = np.round(t[i - 1] + dt, 2)
        #print(Tang)
        #print(Fnet)

        print(Tang)
        print(Fnet)
        print(V[i])
        if t[i] >= duration - dt:
            return A, V, X, t

        while V[i][0] < 400/3.6:
            mode = 3
            i = i + 1
            alpha = np.radians(8)
            acc = 0.5
            W = np.array([-W0 * np.sin(pitch), 0, W0 * np.cos(pitch)])
            Fw, Fh = aerodynamic_force(V[i - 1], alpha, beta, alphah, S)
            T = thrust_force(mode, Fw + Fh, W, acc)
            Fnet = W + Fw + Fh + T
            Fnet = np.array([Fnet[0],0,0])
            A[i] = Fnet / m
            V[i] = V[i - 1] + A[i] * dt
            X[i] = X[i - 1] + V[i] * dt
            t[i] = np.round(t[i - 1] + dt, 2)
            if t[i] >= duration - dt:
                return A, V, X, t


    return A, V, X, t


    # for i,time in enumerate(np.arange(0, duration, dt)):
    #     mode = 0
    #     W = np.array([-W0 * np.sin(pitch), 0, W0 * np.cos(pitch)])
    #     Fw, Fh = aerodynamic_force(V[i-1], alpha, beta, alphah, S)
    #     T = thrust_force(mode, (Fw+Fh), W)
    #     Fnet = W + Fw + Fh + T
    #     A[i] = Fnet / m
    #     V[i] = V[i-1] + A[i] * dt
    #     X[i] = X[i-1] + V[i] * dt
    #     t[i] = time




if __name__ == "__main__":
    A, V, X, t = run_simulation(2000, 0.1)
    print(V.T[2])
    plt.plot(t, X.T[2], color='b')
    Vnew = []
    for i in V:
        Vnew.append(np.linalg.norm(i))
    plt.plot(t, Vnew, color='r')
    plt.show()





# # Needed:
# # Gravity vector, thrust vector, lift vector, drag vector
# # Inertial tensor
# # Propeller torque-thrust ratio
# labda = 0.85
#
# # Propeller thrust
# T_prop = 5000  # [N]
#
# # Propeller direction
# # d = np.array([cos(w_prop*t)],[0],[sin(w_prop*t)])
# d = np.array([1,0,0])
#
# # Propeller position in body frame
# # r = np.array([x_prop + rotor_length*cos(w_prop*t)],[0],[z_prop + rotor_length*sin(w_prop*t)])
# r = np.array([0.5,0,0.5])
#
# # Mass center of wing in body frame
# s = np.array([[0.5,0,0],
#               [-2.5,0,0]])
#
# # Number of propellers = 2
# N_props = 2
#
# # Number of wings = 2
# M_wings = 2
#
# # Air density = 0.01
# rho_Mars = 0.01
#
# # Angle of attack = 0
# alpha0 = 0
#
# # Surface area of each wing
# # State vectors of position, Euler angles, velocity and angular velocity