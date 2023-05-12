

import np as np
from matplotlib import pyplot as plt

from dse.preliminary.AircraftEstimating import Class2Weight
#Payload Range Diagram


def PayloadRange():
    payloadmass=np.arange(400, 689, 1)
    RangeArr = Class2Weight(Wto=3000, N_ult=4.4, AR=28, wingbraced=True, V_cr=154, E_density=333, m_payload=payloadmass, m_solar=0)

    RangeLimity, RangeLimitx = [0, 900], [1000, 1000]
    PayloadLimity, PayloadLimitx = [400, 400], [0, 6000]
    PayRange, axis = plt.subplots()

    axis.plot(RangeArr, payloadmass)
    axis.plot(RangeLimitx, RangeLimity, label='Limit set by minimum Range requirement')
    axis.plot(PayloadLimitx, PayloadLimity, label='Limit set by minimum Payload requirement')
    plt.legend(loc='best')
    plt.xlim(0, 6000)
    plt.ylim(250, 800)
    plt.savefig('Tiltrotor-PayloadRangeDiagram.pdf')
    plt.show()


#Climb Performance
def AircraftClimbPerf(T, W, Cl, Cd):
    gamma_climb = np.arctan(T/W - Cd/Cl)

def RotorClimbPerf(MTOM, R):
    MTOW = MTOM * 0.45
    A = np.pi*R**2
    T_max= 1.1*MTOM*3.721/4
    T_hover = MTOM*3.721/4
    rho=0.01
    V_c = 0*3.28084 #Climb Velocity
    v_1max = np.sqrt(T_max/(2*rho*A)) #m/s
    v_1hover = np.sqrt(T_hover/(2*rho*A)) #m/s
    v_1c = -V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)

    DhpActual = T_max*v_1max - T_hover*v_1hover /745.7
    Dhp=0
    v_1hover *= 3.28084
    while Dhp <= DhpActual:
        V_c+=0.1
        Dhp = MTOW/550 *(V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)-v_1hover)

    print(V_c/3.28084)
    #Dhp = (MTOW*(v_1c+V_c)+4*(Dv/MTOW)*rho/2*(v_1c+V_c)**3*A_M+(DA*Cd)*)



#Calling Previous functions
PayloadRange()
RotorClimbPerf(3000, 14.1)