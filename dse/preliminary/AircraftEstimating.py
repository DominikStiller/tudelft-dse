import numpy as np
import scipy.integrate
from scipy.interpolate import InterpolatedUnivariateSpline
from RotorEngineSizing import RadiusMassElementMomentum
from PerformanceAnalysis import RangeCalc
from cruise_sizing import area_and_thrust

Wcrew = 300 * 2.205 #kg
Wpayload = 100 * 2.205 #kg

we_w0 = 0.55

R = 2000*1000*3.28
C = 1/3600 # kg/hr
V = 150*3.28 #m/s
L2D = 1.5/0.11

wc_w0 = np.exp(-R*C/(V * L2D))

wf_w0 = 1.06*(1-wc_w0)

A=0.72
C=-0.03
w0_init = 400
w0_found = 300
while (w0_init-w0_found)/w0_found > 0.01:
    we_w0 = A * w0_init ** C

    w0_found = (Wcrew + Wpayload) / (1 - wf_w0 - we_w0)
print(w0_found/2.205)


print(wf_w0*w0_found/2.205)

def DragEstimation(lf, hf, Swing, t2c, Vcr, visc_cr, Cl, AR, rho):
    Oswald = 0.9
    #Fuselage
    bf=hf
    Cd_fusS = 0.0031*lf*(bf+hf)

    #Wing
    Cd_wingS = 0.0054*(1+3*t2c*np.cos(0)**2)*Swing #Sweep is 0

    #Engines
    Cd_engineS = 2*0.07*0.08

    #Undercarriage
    r_uc = 1.08 #Main gear retracted in strealined fairings
    #Reynolds correction
    r_Re = 47*(Vcr*lf/visc_cr)**(-0.2)
    #Tailplane
    r_t = 1.24

    Cd_0 = r_Re*r_uc*(r_t*(Cd_fusS+Cd_engineS))/Swing

    Cd = Cd_0+Cl**2/(np.pi*AR*Oswald)
    D = Cd*rho*0.5*Vcr**2*Swing
    print('Drag of the aircraft without wing will be of '+str(D)+'[N]')
    return D



def RangeCalc(Wto, Wtot, R, AR, V_cr, E_density):
    g_mars=3.721
    b = 2 * 1.5 * R
    c = b / AR
    bf = c
    lf = 1.78 + c * 3
    if c < 1.78:
        lf = 1.78 * 4
        bf = 1.78
    hf=bf
    Swing = (b - bf) * c
    T_to = 1.1 * Wto * g_mars
    T_cr = area_and_thrust(0, 1.2, 0.11, Wto, 0.5 * 0.01 * V_cr)[1] + DragEstimation(lf, hf, Swing, 0.12, V_cr,
                                                                                     5.167E-4, 1.2, AR, rho=0.01)
    Power = T_to * np.sqrt(T_to / (2 * 0.01 * np.pi * R * R))
    print(T_cr * np.sqrt(T_cr / (2 * 0.01 * np.pi * R * R)))
    E_TO = Power * (4 / 6)
    E_cr = ((Wto - Wtot - 400) * E_density - E_TO)
    P_cr = T_cr * np.sqrt(T_cr / (2 * 0.01 * np.pi * R * R))
    Endurance = E_cr / P_cr
    Range = Endurance * V_cr * 3.6
    return Range, Endurance


#Weight Prediction:
def Class2Weight(Wto, N_ult, AR, wingbraced, V_cr, E_density):
    g_mars=3.721
    R, T, hp, Power, RotorMass = RadiusMassElementMomentum(Wto, 4, coaxial=True)
    b=2*1.5*R
    c=b/AR
    bf=c
    lf = 1.78+c*3
    if c<1.78:
        lf=1.78*4
        bf=1.78
    Swing = (b-bf)*c
    Stail = Swing*c/(1.5*R)

    hf=bf

    #Crude Estimation
    Wstruc2Wto = 0.447*np.sqrt(N_ult)*(bf*hf*lf/Wto)**0.24


    #Wing Group
    Wwing2Wto = 4.9e-3*b**0.75*(1+np.sqrt(1.905/b))*N_ult**0.55*((b/c)/(Wto/Swing))**0.3
    if wingbraced == True:
        Wwing2Wto*=0.7

    #Tail Group
    Wtail2Wto = 0.64*(N_ult*Stail**2)**0.75

    #Body Group
    lamdaf = lf/hf
    Vdive = 1.25*V_cr
    lt = lf
    S_g = np.pi*hf*lf*(1-2/lamdaf)**(2/3)*(1+1/(lamdaf**2))
    Wf = .23*np.sqrt(Vdive*lt/(bf+hf))*S_g**1.2

    #Control Surfaces
    ksc=0.64 #transport plane with powered controls
    Wsc = 0.768*ksc*Wto**(2/3)

    #Total Weight
    Wtot = Wwing2Wto*Wto+Wtail2Wto+Wf+Wsc+RotorMass
    Range, Endurance = RangeCalc(Wto, Wtot, R, AR, V_cr, E_density)

    print('Very Crude Structural Estimate: '+str(Wstruc2Wto*Wto))
    print('Less crude estimate: \n')
    print('Wing weight: '+str(Wwing2Wto*Wto)+'[kg]')
    print('Tail weight: ' + str(Wtail2Wto) + '[kg]')
    print('Body weight: '+str(Wf)+'[kg]')
    print('Control Surfaces: '+str(Wsc)+'[kg]')
    print('Available weight for batteries: '+str(Wto-Wtot-400)+'[m]')
    print('Available Endurance: '+str(Endurance)+'[h]')
    print('Available Range: '+str(Range)+'[km]')
    print('Flight Radius: '+str(Range/2)+'[km]')

Class2Weight(Wto=3000, N_ult=4.4, AR=28, wingbraced=True, V_cr=154, E_density=333)

