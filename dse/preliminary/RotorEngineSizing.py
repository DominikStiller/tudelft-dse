import numpy as np
import matplotlib.pyplot as plt
#Given Values
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline


def RadiusfromMass(M):
    N_rotor=4
    T_init = 1.1*M * 3.721 / N_rotor
    b=6 #Number of blades per rotor
    a=6 #Lift Curve slope [1/rad]
    R=1
    cd=0.06
    theta_tip = 10 *np.pi/180
    T=0

    def func(x):
        return theta_tip*R*x - x**2 * np.arctan(v1/(omega*x))

    def func2(x):
        return x*c*rho/2*omega**2 * (x**2*a*(theta_tip/(x/R)-np.arctan(v1/(omega*x)))*np.arctan(v1/(omega*x))+cd*x**2)

    while T<T_init:
        R=R+0.1
        c = R / 25  # Chord of each rotor
        A = np.pi * R ** 2  # Area of the rotor
        omega = 200 / R
        v1 = np.sqrt(T_init / (2 * rho * A))
        #T=0.88 * b * rho * omega**2 * a *c/2*(theta_tip * R**3 /2-(v1 * (np.log(v1**2 + (omega*R)**2)-np.log(v1**2))/(2*omega)))
        integral, accuracy = integrate.quad(func, 0, R)
        T = 0.88 * b/2*rho*omega**2 *a*c * integral
        Q, Q_accuracy = integrate.quad(func2, 0, R)
        Q = Q*0.88*b
    sigma = b*c/(np.pi*R)
    Cq = Q/(rho*sigma*A*(omega*R)**2*R)
    Ct = T/(rho*sigma*A*(omega*R)**2)
    FM = np.sqrt(sigma/2) * Ct**(3/2)*Cq
    Hp = T*v1*N_rotor*0.00134
    return R, accuracy, Hp, sigma, Ct, Cq, FM

def RadiusMassElementMomentum(M, N_rotors, coaxial):
    gm=3.721
    T_min = 1.1*M*gm /N_rotors
    N_blades = 6
    b=N_blades
    T=0
    R=1
    while T<T_min:
        R+=0.1
        c=R/20
        theta_tip = 6*np.pi/180
        x0 = 0.1
        V_m = 220
        v_tip = 0.92*V_m
        omega = v_tip/R
        rho = 0.01
        n_elements = 10
        a0 = 6
        alpha0 = 2*np.pi/180
        A = np.pi*R**2
        r2R = np.arange(1, n_elements+1)/n_elements
        c2R = c/R
        M_local = (r2R)*(omega*R/V_m)
        a = a0/(np.sqrt(1-M_local**2))
        Dtheta = theta_tip/r2R
        theta0 = -min(Dtheta-np.ones(n_elements,)*alpha0)+theta_tip
        theta = theta0+Dtheta-alpha0
        v12Omegar = a*b*c2R/(16*np.pi*r2R)*(-1+np.sqrt(1+(32*np.pi*theta*r2R)/(a*b*c2R)))
        alpha = theta-np.arctan(v12Omegar)
        alpha = 6*np.pi/180 * np.ones(n_elements)
        cl = a*alpha
        cd = 0.025
        DctDr2R = b*r2R**2*c2R*cl/(2*np.pi)

        funCT = InterpolatedUnivariateSpline(r2R, DctDr2R)
        Ct_lossless = funCT.integral(0, 1)
        if Ct_lossless<0.006:
            B = 1-0.06/b
        else: B = 1-np.sqrt(2.27*Ct_lossless-0.01)/b

        Ct = Ct_lossless - funCT.integral(B, 1)

        Dcq0Dr2R = b*r2R**3*c2R*cd/(2*np.pi)
        funCQ0 = InterpolatedUnivariateSpline(r2R, Dcq0Dr2R)
        CQ_profile = funCQ0.integral(x0, 1)
        DcqiDr2R = b*r2R**3*c2R*cl*v12Omegar/(2*np.pi)
        funCQi = InterpolatedUnivariateSpline(r2R, DcqiDr2R)
        CQ_induced = funCQi.integral(B, 1)

        DCQ_I = 0.01*CQ_induced

        DL = Ct*rho*(omega*R)**2

        Ct2sigma = Ct/(b*c/(np.pi*R))

        Cq = (CQ_profile+CQ_induced+DCQ_I)*.93

        T=rho*A*(omega*R)**2*Ct
        if coaxial == True:
            T*=0.88
        power = rho*A*(omega*R)**3*Cq

    #Rotor Weight:
    def NACA0012airfoilfunction(x):
        return 5 * 0.12 * (0.2969* (x**0.5) - 0.1260 * x - 0.3516 * (x ** 2) + 0.2843 * (x ** 3) - 0.1036 * (x ** 4))


    Area = integrate.quad(NACA0012airfoilfunction, 0, 1)[0] * c
    fillfactor=0.3
    Rotor_mass = 1500*R*Area*fillfactor

    print('Given the mass of '+str(M)+'kg, the following applied: \n')
    print('Rotor Radius: '+str(R)+'[m]')
    print('Rotor Power: '+str(power)+'[W]')
    print('Total Aircraft power: '+str(power*N_rotors)+'[W]')
    print('Weight per rotor with '+str(b)+' blades is '+ str(Rotor_mass)+'[kg]')

    return R, T, power/745.7, power

if __name__ == '__main__':
    gm = 3.721
    MTOW = 3000 * gm
    rho = 0.01  # kg/m3

    # Blade Element Method with ideal twist
    N_rotor = 4

RadiusMassElementMomentum(3000, 4, coaxial=True)




''''
hp = T*v1*0.01315
DL = T/A

plt.scatter(R, T)
plt.plot(R, T_init*np.ones(np.shape(R)))
plt.show()

#Using above values:

R1=21.5#[m]
A1=np.pi* R1**2
omega1 = 200/R1
c1=R1/25
sigma = 6*c1/(np.pi*R1)
Ct = T_init/(rho*A1*(omega1*R1)**2)

T=np.empty(23,)
hp=np.empty(23,)
DL=np.empty(23,)
for R in range(2, 25):
    c = R / 25
    A = np.pi * R ** 2
    omega = 200 / R
    v1 = np.sqrt(T_init / (2 * rho * A))
    theta_tip = 10 * np.pi / 180
    T[R-2] = b * rho * omega ** 2 * a * c / 2 * (
                theta_tip * R ** 3 / 2 - (v1 * (np.log(v1 ** 2 + (omega * R) ** 2) - np.log(v1 ** 2)) / (2 * omega)))
    hp[R-2] = T[R-2] * v1 * 0.01315
    DL[R-2] = T[R-2] / A

R=np.arange(2, 25, 1)

plt.scatter(R, T)
plt.plot(R, T_init*np.ones(np.shape(R)))
plt.show()

plt.scatter(R, hp)
plt.show()

plt.scatter(R, DL)
plt.show()
'''