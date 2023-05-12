import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from cruise_sizing import max_tipSpeed
from power_sizing import power
from constants import const


def RadiusMassElementMomentum(M, N_rotors, N_blades, coaxial, V_tip, print_results=False):
    # Rename constants
    gm = const['gravityMars']
    rho = const['airDensity']
    V_m = const['soundSpeed']

    b = N_blades
    v_tip = V_tip
    T_min = 1.1 * M * gm /N_rotors

    # Initial guesses
    T = 0
    R = 1
    while T < T_min:
        R += 0.1  # Increase the radius
        c = R/20  # Update the chord

        theta_tip = np.radians(6)  # Assumed angle
        x0 = c/R  # Distance from the blade's base to the centre of rotation
        omega = v_tip/R  # Angular velocity

        n_elements = 10
        a0 = 6
        alpha0 = np.radians(2)
        A = np.pi*R**2

        r2R = np.arange(1, n_elements+1)/n_elements
        c2R = c/R
        M_local = (r2R)*(omega*R/V_m)  # Local Mach number (omega*R = V_tip which is constant)
        a = a0/(np.sqrt(1-M_local**2))
        Dtheta = theta_tip/r2R
        theta0 = -min(Dtheta-np.ones(n_elements,)*alpha0)+theta_tip
        theta = theta0+Dtheta-alpha0
        v12Omegar = a*b*c2R/(16*np.pi*r2R)*(-1+np.sqrt(1+(32*np.pi*theta*r2R)/(a*b*c2R)))

        alpha = 6*np.pi/180 * np.ones(n_elements)

        #S1223
        cl = 1.6
        cd = 0.05
        DctDr2R = b*r2R**2*c2R*cl/(2*np.pi)

        funCT = InterpolatedUnivariateSpline(r2R, DctDr2R)
        Ct_lossless = funCT.integral(x0, 1)
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

        T = rho*A*(omega*R)**2*Ct
        if coaxial:
            T *= 0.88

        pow = power(T, R)

    #Rotor Weight:
    def NACA0012airfoilfunction(x):
        return 5 * 0.12 * (0.2969 * (x**0.5) - 0.1260 * x - 0.3516 * (x ** 2) + 0.2843 * (x ** 3) - 0.1036 * (x ** 4))


    Area = integrate.quad(NACA0012airfoilfunction, 0, 1)[0] * c
    fillfactor = 0.3
    Rotor_mass = 1500*R*Area*fillfactor

    if print_results:
        print('Given the mass of '+str(M)+'kg, the following applied: \n')
        print('Rotor Radius: '+str(R)+'[m]')
        print('Rotor Power: ' + str(pow) + '[W]')
        print('Total Aircraft power: ' + str(pow * N_rotors) + '[W]')
        print('Weight per rotor with '+str(b)+' blades is '+ str(Rotor_mass)+'[kg]')
        print('Total Rotor mass for entire aircraft: '+str(Rotor_mass*N_rotors))

    return R, T, pow / 745.7, pow * N_rotors, Rotor_mass * N_rotors


if __name__ == '__main__':
    gm = const['gravityMars']
    rho = const['airDensity']  # kg/m3
    MTOW = 3000 * gm

    # Blade Element Method with ideal twist
    N_rotor = 4

    RadiusMassElementMomentum(3000, 4, 6, coaxial=True,V_tip= 168)

    cruiseSpeed = np.linspace(112, 154)
    a, b = max_tipSpeed(cruiseSpeed)
    V_tip = a * cruiseSpeed + b

    R = np.empty(np.shape(V_tip))
    for i in range(len(V_tip)):
        R[i] = RadiusMassElementMomentum(3000, 4, 6, coaxial=True, V_tip=V_tip[i])[0]

    def oneOverX(x, a, b, c, d):
        return a/(b*x+c) + d

    popt2, pcov2, infoDict, msg, ier = curve_fit(oneOverX, V_tip, R, full_output=True)
    print(f'Rotor Radius =  {popt2[0]}/({popt2[1]}x + {popt2[2]}) + {popt2[3]}')
    print(pcov2)
    plt.scatter(V_tip, R, alpha=0.5)
    plt.plot(V_tip, oneOverX(V_tip, *popt2), 'r--')
    plt.xlabel('Tip speed [m/s]')
    plt.ylabel('Rotor radius [m]')
    plt.show()


    plt.plot(cruiseSpeed, oneOverX(a*cruiseSpeed+b, *popt2), label=f'{np.round(popt2[0], 2)}/'
                                                                   f'({np.round(popt2[1]*a, 2)}cruiseSpeed + '
                                                                   f'{np.round(popt2[1]*b + popt2[2], 2)}) + '
                                                                   f'{np.round(popt2[3], 2)})')
    plt.xlabel('Cruise speed [m/s]')
    plt.ylabel('Rotor radius [m]')
    plt.legend()
    plt.show()
    print(f'rotorRadius = {np.round(popt2[0], 2)}/({np.round(popt2[1]*a, 2)}cruiseSpeed + '
          f'{np.round(popt2[1]*b + popt2[2], 2)}) + {np.round(popt2[3], 2)})')

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