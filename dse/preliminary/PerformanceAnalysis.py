

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

def RotorClimbPerf():
    V_c = 3*3.28084 #Climb Velocity
    v_1hover = 14*3.28084 #m/s
    MTOW = 3000 * 0.45
    v_1c = -V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)

    Dhp = MTOW/550 *(V_c/2 + np.sqrt((V_c/2)**2+v_1hover**2)-v_1hover)
    print(Dhp)
    #Dhp = (MTOW*(v_1c+V_c)+4*(Dv/MTOW)*rho/2*(v_1c+V_c)**3*A_M+(DA*Cd)*)



#Calling Previous functions
PayloadRange()