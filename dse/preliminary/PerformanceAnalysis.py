

#Payload Range Diagram

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

#Climb Performance
import numpy as np

V_c = 3*3.28084 #Climb Velocity
v_1hover = 14*3.28084 #m/s
MTOW = 3000 * 0.45
v_1c = -V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)

Dhp = MTOW/550 *(V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)-v_1hover)
print(Dhp)
#Dhp = (MTOW*(v_1c+V_c)+4*(Dv/MTOW)*rho/2*(v_1c+V_c)**3*A_M+(DA*Cd)*)

#Low Rate Descent
Vd = 10#descent velocity
v_1d = Vd/2 + np.sqrt((Vd/2)**2+v_1hover**2)
#Vortex Ring

Vd = 2*v_1hover
#windmill brake
v_1wb = Vd/2 - np.sqrt((Vd/2)**2-v_1hover**2)

#Take-off and Landing
