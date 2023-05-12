import numpy as np
from matplotlib import pyplot as plt

from dse.preliminary.AircraftEstimating import Class2Weight
#Payload Range Diagram


def PayloadRange(R, mass_rotors, Mass_design, N_ult, AR, wingbraced, V_cr, E_density, P_density_TO, E_density_TO, Mass_solar, minpayload, minRange):
    payloadmass=np.arange(minpayload, 700, 1)
    RangeArr = Class2Weight(R, mass_rotors, Mass_design, N_ult, AR, wingbraced, V_cr, E_density, P_density_TO, E_density_TO, payloadmass, Mass_solar)[0]

    RangeLimity, RangeLimitx = [0, 900], [minRange, minRange]
    PayloadLimity, PayloadLimitx = [400, 400], [0, 20000]
    PayRange, axis = plt.subplots()

    axis.plot(RangeArr, payloadmass)
    axis.plot(RangeLimitx, RangeLimity, label='Limit set by minimum Range requirement')
    axis.plot(PayloadLimitx, PayloadLimity, label='Limit set by minimum Payload requirement')
    plt.legend(loc='best')
    plt.ylim(250, 800)
    plt.savefig('Tiltrotor-PayloadRangeDiagram.pdf')
    plt.show()


#Climb Performance
def AircraftClimbPerf(m_bat_cr, P_dens_cr, Cl, Cd, M_design, R, V_cr):
    W=M_design*3.71
    P_climb_ac = m_bat_cr * P_dens_cr
    #Climb angle in cruise
    T_max = ( P_climb_ac/(np.sqrt(1 / (2 * 0.01 * np.pi * R * R))) )**(2/3)
    gamma = np.arctan(T_max/W -Cd/Cl)
    ROC_cr = np.sin(gamma)*V_cr
    return ROC_cr

def RotorClimbPerf(MTOM, R):
    MTOW = MTOM * 0.45
    A = np.pi*R**2
    T_max= 1.1*MTOM*3.71/4
    T_hover = MTOM*3.71/4
    rho=0.01
    V_c = 0*3.28084 #Climb Velocity
    v_1max = np.sqrt(T_max/(2*rho*A)) #m/s
    v_1hover = np.sqrt(T_hover/(2*rho*A)) #m/s
    v_1c = -V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)

    DhpActual = (T_max*v_1max - T_hover*v_1hover) /745.7
    Dhp=0
    v_1hover *= 3.28084
    while Dhp <= DhpActual:
        V_c+=0.1
        Dhp = MTOW/550 *(V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)-v_1hover)
    ROC_vert = V_c

    return ROC_vert
