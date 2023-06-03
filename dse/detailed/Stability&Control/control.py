import matplotlib.pyplot as plt
import numpy as np
import scipy

I = [20000, 20000, 20000] #moment of inertia around x, y, z axis
rho = 0.01
clalpha = 2 * np.pi
cl0 = 1
clalphah = 2 * np.pi
clcd = 20
clcdh = 20
S = [100, 10]
m = 2700
gmars = 3.71
Vclimb = 2
Pmax = 34000
Rrotor = 10.4
Vhover = 19.52

#simulation
duration = 20 #duration of simulation in seconds
dt = 0.1 #time increment in seconds
Tangvel = np.radians(4.5*dt) #3 deg/s ang vel of thrust

def definegeometry(mode, Tang):
    '''

    Returns: the position of each force with respect to the aircraft body reference frame

    '''
    if mode == 0:
        # Thrust is straight up
        RThrustLeft = np.array([1,-17,-2])
        RThrustRight = np.array([1,17,-2])
    elif mode == 1:
        # Thrust is straight up
        RThrustLeft = np.array([1,-17,-2])
        RThrustRight = np.array([1,17,-2])
    elif mode == 2:
        # Thrust is turning with Tang
        Tang = np.radians(Tang)
        RThrustLeft = np.array([1+np.cos(abs(Tang)),-17,-1-np.sin(abs(Tang))])
        RThrustRight = np.array([1+np.cos(abs(Tang)),17,-1-np.sin(abs(Tang))])
    elif mode == 3:
        # Thrust is straight ahead
        RThrustLeft = np.array([2,-17,-1])
        RThrustRight = np.array([2,17,-1])
    Rwingleft = np.array([1, -10, -1])
    Rwingright = np.array([1, 10, -1])
    Rstabilizer = np.array([-10, 0, -2])
    R = np.array([Rwingleft, Rwingright, Rstabilizer, RThrustLeft, RThrustRight])
    return R

def aero_to_body(alpha, beta):
    """
    Transformation matrix of aerodynamic reference frame to body reference frame
    :param alpha: Angle of attack, alpha, [deg]
    :param beta: Sideslip angle, beta, [deg]
    :return: Transformation matrix
    """
    a = np.radians(alpha)
    b = np.radians(beta)
    T = np.array([[np.cos(b) * np.cos(a), np.sin(b), np.cos(b) * np.sin(a)],
         [-np.sin(b) * np.cos(a), np.cos(b), -np.sin(b) * np.sin(a)],
         [-np.sin(a), 0, np.cos(a)]])
    T = np.linalg.inv(T)
    return T

def body_to_inertial(theta, gamma, psi):
    """
    Transformation matrix of body reference frame to inertial reference frame
    :param theta: Roll angle, theta, [deg]
    :param gamma: Flight path angle, gamma, [deg]
    :param psi: Heading angle, psi, [deg]
    :return: Transformation matrix
    """
    t = np.radians(theta)
    g = np.radians(gamma)
    p = np.radians(psi)
    T = np.array([[np.cos(t) * np.cos(p), np.cos(t) * np.sin(p), -np.sin(t)],
                  [np.sin(g) * np.sin(t) * np.cos(p) - np.cos(g) * np.sin(p), np.sin(g) * np.sin(t) * np.sin(p) + np.cos(g) * np.cos(p), np.sin(g) * np.sin(t)],
                  [np.cos(g) * np.sin(t) * np.cos(g) + np.sin(g) * np.sin(p), np.cos(g) * np.sin(t) * np.sin(p) - np.sin(g) * np.cos(p), np.cos(g) * np.cos(t)]])
    T = np.linalg.inv(T)
    return T

def thrust_force(mode, Faero, W, ang, V):
    if mode == 0:
        Tmax = 12000
        T = np.array([0, 0, -Tmax])
        return T
    elif mode == 1:
        T = -1 * (W + Faero)
        return T
    elif mode == 2:
        # T = -1 * W - np.array([0, 0, Faero[2]])
        # newang = ang + Tangvel
        # T[0] = T[2] / np.tan(newang)
        # print(T[0])

        # Tmax = min((Pmax*4 / (0.1047198*183*5.2 - np.sin(ang) * V[0])), 12000)
        # V1 = np.sqrt(-(V[0]**2)/2 + np.sqrt(((V[0]**2)/2)**2 + Vhover**4))
        # Trotor_1 = rho * np.pi * Rrotor ** 2 * np.sqrt(V[0] ** 2 + V1 ** 2) * 2 * V1
        # Tmax_1 = Trotor_1 * 4
                # V1 = (omega * Rrotor)/mu * (CT / 2)
        # Trotor = rho * np.pi * Rrotor ** 2 * np.sqrt(V[0]**2 + V1**2) * 2 * V1
        if V[0] == 0:
            Tmax = 12000
        else:
            omega = 19.23
            CT = 0.02239
            mu = V[0] / (omega * Rrotor) * np.cos(ang)
            Cd = 0.05
            Crotor = 1 / 20 * Rrotor
            Trotor_2 = np.sqrt((Pmax + Faero[0] / 2 * V[0] - (rho * 6 * Rrotor * Crotor * ((omega * Rrotor) ** 3) * Cd * (
                        1 + 3 * mu ** 2)) / 8) * rho * np.pi * Rrotor ** 2 * V[0])
            Tmax = Trotor_2 * 4
        # if V[0] != 0 and Tmax_2 <= Tmax_1:
        #     Tmax = Tmax_2
        # elif V[0] != 0 and Tmax_2 > Tmax_1:
        #     Tmax = Tmax_1
        print(f"AAAAAAAHHHHHHHH {Tmax}, {ang*180/np.pi}, {V[0]}, {Faero[0]}")
        # print(V[0])
        # print(Tmax)
        T = np.array([np.cos(ang) * Tmax, 0, -np.sin(ang)*Tmax])
        return T
        # if abs(Tmax) > np.linalg.norm(T):
        #     print(f'Thrust less than max thrust')
        #     return T, newang
        # else:
        #     print(f'Thrust higher than max thrust')
        #     print(T[0])
        #     T[0] = np.sqrt(Tmax**2 - T[2]**2)
        #     print(T[0])
        #     return T, ang
    elif mode == 3:
        if V[0] > 400/3.6:
            Tcurrent = Faero[0] - m * ang
            T = np.array([Tcurrent, 0, 0])
            return T
        else:
            Tcurrent = Faero[0] + m * ang
            T = np.array([Tcurrent, 0, 0])
            return T

def aerodynamic_force(V, alpha, beta, alphah, S, imain):
    alpha = np.radians(alpha)
    beta = np.radians(beta)
    alphah = np.radians(alphah)
    imain = np.radians(imain)
    CLw = clalpha * alpha + cl0
    CDw = CLw / clcd

    CLh = clalphah * alphah
    CDh = CLh / clcdh

    Cdwz = 0.2
    Cdhz = 0.2

    pitch = alpha - imain

    q = 0.5 * rho

    Fwz = - CLw * q * S[0] * V[0] ** 2 - (Cdwz * q * S[0] * V[2] ** 2)
    Fhz = - CLh * q * S[1] * V[0] ** 2 - (Cdhz * q * S[1] * V[2] ** 2)

    Fwx = - CDw * q * S[0] * V[0] ** 2
    Fhx = - CDh * q * S[1] * V[0] ** 2

    Fw = np.array([Fwx, 0, Fwz]).T
    Fh = np.array([Fhx, 0, Fhz]).T

    Fw = np.matmul(aero_to_body(pitch, beta), Fw).T
    Fh = np.matmul(aero_to_body(pitch, beta), Fh).T

    return Fw/2, Fw/2, Fh

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
        acceleration of X, Y, Z axis [m/s^2]
    """
    sumX = []
    sumY = []
    sumZ = []
    for i in range(len(F)):
        sumX.append(F[i][0])
        sumY.append(F[i][1])
        sumZ.append(F[i][2])

    ax = np.array(sumX)/m
    ay = np.array(sumY)/m
    az = np.array(sumZ)/m
    acc = [ax, ay, az]
    return acc

def Moment(F,R):
    """
    Calculates moment based on N forces and their position
    :param F: Force array of dimension (3,N)
    :param R: Position array of dimension (3,N)
    :return: Array of dimension (3,1) with total moment around x, y and z
    """
    M = np.zeros((len(R), 3))
    for i in range(len(R)):
        M[i] = np.cross(R[i], F[i])
    M = np.array([sum(M.T[0]), sum(M.T[1]), sum(M.T[2])])
    return M

def run_simulation(duration, dt):
    pitch = 0
    alpha = 0    # will be a function of pitch and the velocity in body
    alphah = alpha # will be a function of alpha
    beta = 0
    W0 = m * gmars
    t = np.zeros(int(duration / dt))
    Tang = np.zeros(int(duration / dt))
    V = np.zeros((int(duration / dt), 3))
    X = np.zeros((int(duration / dt), 3))
    A = np.zeros((int(duration / dt), 3))
    T = np.zeros((int(duration/dt), 3))
    Vi = np.array([0, 0, 0])
    Xi = np.array([0, 0, 0])
    Ai = np.array([0, 0, 0])
    Ti = np.array([0, 0, 0])
    ti = 0
    V[0] = Vi
    X[0] = Xi
    A[0] = Ai
    T[0] = Ti
    t[0] = ti
    i = 0
    # initial acceleration phase
    print(f"Start accelerate to climb speed")
    while V[i][2] > -Vclimb:
        i = i + 1
        mode = 0
        W = np.array([-W0 * np.sin(pitch), 0, W0 * np.cos(pitch)])
        Fwl, Fwr, Fh = aerodynamic_force(V[i-1], alpha, beta, alphah, S, 0)
        Fw = Fwl + Fwr
        T[i] = thrust_force(mode, (Fw + Fh), W, 0, V[i-1])
        R = definegeometry(mode,90)
        # M = Moment(np.array([Fwl, Fwr, Fh, T/2, T/2]), R)
        Fnet = W + Fw + Fh + T[i]
        # print(W)
        # print(Fw)
        # print(Fh)
        # print(T)
        A[i] = Fnet / m
        V[i] = V[i - 1] + A[i] * dt
        # print(V[i][0])
        X[i] = X[i - 1] + V[i] * dt
        t[i] = np.round(t[i - 1] + dt, 2)
        if t[i] >= duration - dt:
            return A, V, X, t, Tang, T

    # climb until 500m phase
    print(f'Climb until 500m')
    while X[i][2] > -500:
        i = i + 1
        mode = 1
        W = np.array([-W0 * np.sin(pitch), 0, W0 * np.cos(pitch)])
        Fwl, Fwr, Fh = aerodynamic_force(V[i - 1], alpha, beta, alphah, S, 0)
        Fw = Fwl + Fwr
        T[i] = thrust_force(mode, Fw+Fh, W, 0, V[i-1])
        Fnet = W + Fw + Fh + T[i]
        R = definegeometry(mode, 90)
        # M = Moment(np.array([Fwl, Fwr, Fh, T/2, T/2]), R)
        A[i] = Fnet / m
        V[i] = V[i - 1] + A[i] * dt
        X[i] = X[i - 1] + V[i] * dt
        t[i] = np.round(t[i - 1] + dt, 2)
        if t[i] >= duration - dt:
            return A, V, X, t, Tang, T[i]

    V[i][2] = 0
    Tang[0:i+1] = np.radians(90)

    print(f"Start with transition")
    currentcounter = i
    while np.round(Tang[i],4) > 0.006:
        i = i + 1
        mode = 2
        imain = 4
        pitch = 0
        alpha = pitch + imain
        # Tang[i] = np.pi/2 - (Tangvel * dt) * (i - currentcounter)
        W = np.array([-W0 * np.sin(pitch), 0, W0 * np.cos(pitch)])
        # print(V[i-1])
        Fwl, Fwr, Fh = aerodynamic_force(V[i - 1], alpha, beta, alphah, S, imain)
        # print(Fwl)
        # print(Fwr)
        Fw = Fwl + Fwr
        T[i] = thrust_force(mode, Fw + Fh, W, Tang[i-1], V[i-1])
        # print(T)
        # print(T)
        # print(V[i-1])
        # print(Fw)
        R = definegeometry(mode,Tang[i])
        # M = Moment(np.array([Fwl,Fwr,Fh,T/2,T/2]),R)
        Fnet = W + Fw + Fh + T[i]
        # print(W)
        # print(Fw)
        # print(Fh)
        # print(T)
        # print(f'Fnet1 = {Fnet}')
        if Fnet[2] <= 0:
            Tang[i] = Tang[i-1] - Tangvel*dt
            T[i][2] = -(W[2] + (Fw[2] + Fh[2]))
            Fnet = W + Fw + Fh + T[i]
        elif Fnet[2] > 0:
            Tang[i] = Tang[i-1]
        print(f'Fnet2 = {Fnet}')
        A[i] = Fnet / m
        V[i] = V[i - 1] + A[i] * dt
        X[i] = X[i - 1] + V[i] * dt
        t[i] = np.round(t[i - 1] + dt, 2)

        if t[i] >= duration - dt:
            return A, V, X, t, Tang, T
    print(f"Start converging to cruise speed")
    # while np.round(V[i][0],2) != np.round(400/3.6,2):
    #     mode = 3
    #     i = i + 1
    #     imain = 4
    #     pitch = 0
    #     alpha = pitch + imain
    #     acc = 0.5
    #     W = np.array([-W0 * np.sin(pitch), 0, W0 * np.cos(pitch)])
    #     Fwl, Fwr, Fh = aerodynamic_force(V[i - 1], alpha, beta, alphah, S, imain)
    #     Fw = Fwl + Fwr
    #     T = thrust_force(mode, Fw + Fh, W, acc, V[i-1])
    #     R = definegeometry(mode,0)
    #     M = Moment(np.array([Fwl,Fwr,Fh,T/2,T/2]),R)
    #     Fnet = W + Fw + Fh + T
    #     Fnet = np.array([Fnet[0],0,0])
    #     A[i] = Fnet / m
    #     V[i] = V[i - 1] + A[i] * dt
    #     X[i] = X[i - 1] + V[i] * dt
    #     t[i] = np.round(t[i - 1] + dt, 2)
    #     if t[i] >= duration - dt:
    #         return A, V, X, t, Tang

    return A[:i], V[:i], X[:i], t[:i], Tang[:i], T[:i]

if __name__ == "__main__":
    A, V, X, t, Tang, T = run_simulation(1000, 0.1)
    plt.figure(1)
    fig, ax1 = plt.subplots()
    color = "tab:red"
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Height [m]", color=color)
    ax1.plot(t, -X.T[2], color=color)
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()

    color = "tab:blue"
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Velocity [m/s]", color=color)
    ax2.plot(t, -V.T[2], color=color)
    ax2.tick_params(axis="y", labelcolor=color)

    fig.tight_layout()
    plt.show()

    Tmax = []
    for i in T:
        Tmax.append(np.linalg.norm(i))
    plt.figure(2)
    plt.plot(V.T[0][int(250/dt):],Tmax[int(250/dt):],color="r")

    plt.show()
    #
    #
    # plt.plot(t, X.T[2], color='b', label='Z Position')
    # plt.plot(t, V.T[2], color='r', label='Z Velocity')
    # plt.legend()
    # plt.figure(2)
    # plt.plot(t,V.T[0], color='g', label='X Velocity')
    # plt.legend()
    # plt.show()


####### scipy.spatial.transform.Rotation ######