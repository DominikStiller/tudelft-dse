import numpy as np
import matplotlib.pyplot as plt
#Given Values
from scipy import integrate

gm=3.721
MTOW = 3000*gm
rho = 0.01 #kg/m3


def Sizing(R):
    Rotor_radius=R
    N_blades = 6
    b = N_blades
    Rotor_chord = 1 / 25 * Rotor_radius
    c = Rotor_chord
    V_sound = 220  # m/s
    omega = 0.8 * V_sound / Rotor_radius
    A = np.pi * Rotor_radius ** 2
    n_elements = 15

    sigma = b * c / (np.pi * Rotor_radius)
    theta_tip = 7 *np.pi/180
    theta_slope = -10*np.pi/180

    theta0 = 3/2 * theta_tip - 3/4 * theta_slope
    alphaL0 = 0
    rtoR=np.empty(n_elements)
    M_local=np.empty(n_elements)
    a=np.empty(n_elements)
    theta=np.empty(n_elements)
    inflowangle=np.empty(n_elements)
    aoa = np.empty(n_elements)
    dCrdrtoR = np.empty(n_elements)
    dCq0drtoR = np.empty(n_elements)
    dCqidrtoR = np.empty(n_elements)

    for i in range(0,n_elements):
        rtoR[i] = (i+1)/n_elements
        nondim_chord = Rotor_chord/Rotor_radius
        M_local[i] = rtoR[i]*((omega*Rotor_radius)/V_sound)
        a[i] = 5.73/(np.sqrt(1-M_local[i])) #Very Crude estimation of lift curve as function of mach

        theta[i] = theta0+theta_slope*rtoR[i]-alphaL0 #rad

        inflowangle[i] = (a[i]*sigma)/(16*rtoR[i]) * (-1 + np.sqrt(1+(32*theta[i]*rtoR[i]/(a[i]*sigma))))
        aoa[i] = (theta[i] - np.arctan(inflowangle[i]))

        cl = aoa[i]*a[i]
        cd = 0.02 - 0.00733333*aoa[i] + 0.00133333 *(aoa[i]**2)

        dCrdrtoR[i] = (rtoR[i]**2*sigma*cl)/2
        dCq0drtoR[i] = b*(rtoR[i]**3)*(nondim_chord)*cd/(2*np.pi)
        dCqidrtoR[i] = b * (rtoR[i] ** 3) * (nondim_chord) * cl * inflowangle[i] / (2 * np.pi)
    print('inflowangle: '+str(theta*180/np.pi))

    derivative = np.polyfit(rtoR, dCrdrtoR, 9)
    derivativeCq0 = np.polyfit(rtoR, dCq0drtoR, 9)
    derivativeCqi = np.polyfit(rtoR, dCqidrtoR, 9)

    def func(x):
        a,b,c,d,e,f,g,h,i,j = derivative
        return a*x + b*x**2 + c*x**3 + d*x**4 + e*x**5 + f*x**6 + g*x**7 + h*x**8 + i*x**9 + j*x**10
    def funcCq0(x):
        a,b,c,d,e,f,g,h,i,j = derivativeCq0
        return a*x + b*x**2 + c*x**3 + d*x**4 + e*x**5 + f*x**6 + g*x**7 + h*x**8 + i*x**9 + j*x**10
    def funcCqi(x):
        a,b,c,d,e,f,g,h,i,j = derivativeCqi
        return a*x + b*x**2 + c*x**3 + d*x**4 + e*x**5 + f*x**6 + g*x**7 + h*x**8 + i*x**9 + j*x**10


    Ct_lossless = np.trapz(dCrdrtoR)
    Ct_lossless = integrate.quad(func, 0, 1)[0]

    Cq0 = integrate.quad(funcCq0, 0, 1)[0]

    if Ct_lossless>0.006:
        B=1-(np.sqrt(2.27*Ct_lossless-0.01)/b)
    else:
        B = 1-0.06/b

    B = 1-np.sqrt(2*Ct_lossless)/b

    Ct = Ct_lossless - integrate.quad(func, B, 1)[0]

    print(Ct/sigma)
    Cqi = integrate.quad(funcCqi, B, 1)[0]
    Cq = Cqi+Cq0
    DiskLoad = Ct*rho*(omega*Rotor_radius)**2
    Thrust = rho*A*(omega*Rotor_radius)**2*Ct
    Horsepower = rho*A*(omega*Rotor_radius)**3*Cq/550
    print('Disk Loading is: '+str(DiskLoad))
    print("Thrust is: "+ str(Thrust))
    print("Horsepower is: "+ str(Horsepower))

    print(Rotor_chord)
    print(omega)

    return Thrust

Radius = np.arange(2, 25, 1)
j=0
Thrust = np.empty(len(Radius))
for Rad in Radius:
    Thrust[j] = Sizing(Rad)

    if Thrust[j]>= 6000:
        print(Thrust[j])
        break
    j = j + 1
plt.plot(Radius, Thrust)
plt.show()


#Blade Element Method with ideal twist
rho=0.01
T_init = 3000*3.721 /2
b=6
a=6
R=np.arange(3, 25, 0.5)
c=R/25
A=np.pi* R**2
omega = 200/R
v1 = np.sqrt(T_init/(2*rho*A))
theta_tip = 10 *np.pi/180
T=b * rho * omega**2 * a *c/2*(theta_tip * R**3 /2-(v1 * (np.log(v1**2 + (omega*R)**2)-np.log(v1**2))/(2*omega)))
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
print(sigma)
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
    #print(v1)
    print(omega*9.5492968)
R=np.arange(2, 25, 1)

plt.scatter(R, T)
plt.plot(R, T_init*np.ones(np.shape(R)))
plt.show()

plt.scatter(R, hp)
plt.show()

plt.scatter(R, DL)
plt.show()