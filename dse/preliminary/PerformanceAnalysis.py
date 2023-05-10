

#Payload Range Diagram
#GW_lndng = GW_to - (UsedFuel-WUTO)
#Payload = GW_to-GW_min - (UsedFuel+ WUTO + Reserves + AuxFuelTank)

#Range = integrate.quad(funcSR, 0, 1)[0]

#Climb Performance
import numpy as np

V_c = 10*3.28084 #Climb Velocity
v_1hover = 14*3.28084 #m/s
MTOW = 3000 * 0.45
v_1c = -V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)

#Dhp = MTOW/550 *(V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)-v_1hover)

Dhp = (MTOW*(v_1c+V_c)+4*(Dv/MTOW)*rho/2*(v_1c+V_c)**3*A_M+(DA*Cd)*)

#Low Rate Descent
Vd = 10#descent velocity
v_1d = Vd/2 + np.sqrt((Vd/2)**2+v_1hover**2)
#Vortex Ring

Vd = 2*v_1hover
#windmill brake
v_1wb = Vd/2 - np.sqrt((Vd/2)**2-v_1hover**2)

#Take-off and Landing
