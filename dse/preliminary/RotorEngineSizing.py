import numpy as np
import matplotlib.pyplot as plt
#Given Values
from scipy import integrate

gm=3.721
MTOW = 3000*gm
rho = 0.01 #kg/m3

#Blade Element Method with ideal twist
N_rotor = 4
def RadiusfromMass(M):
    T_init = M * 3.721 / N_rotor
    b=6 #Number of blades per rotor
    a=6 #Lift Curve slope [1/rad]
    R=1
    c=R/25 #Chord of each rotor
    A=np.pi* R**2 #Area of the rotor
    omega = 200/R
    v1 = np.sqrt(T_init/(2*rho*A))
    theta_tip = 10 *np.pi/180
    T=0

    def func(x):
        return theta_tip*R*x - x**2 * np.arctan(v1/(omega*x))

    while T<T_init:
        R=R+0.1
        #T=0.88 * b * rho * omega**2 * a *c/2*(theta_tip * R**3 /2-(v1 * (np.log(v1**2 + (omega*R)**2)-np.log(v1**2))/(2*omega)))
        integral, accuracy = integrate.quad(func, 0, R)
        T = 0.88 * b*rho*omega**2 *a*c * integral
    return R, accuracy

print(RadiusfromMass(3000))
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