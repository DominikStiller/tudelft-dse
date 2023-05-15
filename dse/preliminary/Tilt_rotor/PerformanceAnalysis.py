import numpy as np
from constants import const
from matplotlib import pyplot as plt
from AircraftEstimating import Class2Weight

#Payload Range Diagram


def PayloadRange(R, mass_rotors, Mass_design, N_ult, AR, wingbraced, V_cr, E_density, P_density_TO, E_density_TO, Mass_solar, minpayload, minRange, availableSpace):
    payloadmass = np.arange(minpayload, minpayload+availableSpace)
    RangeArr = Class2Weight(R, mass_rotors, Mass_design, N_ult, AR, wingbraced, V_cr, E_density, P_density_TO, E_density_TO, payloadmass, Mass_solar, print_results=False)[0]

    RangeLimity, RangeLimitx = [0, max(payloadmass)], [minRange, minRange]
    PayloadLimity, PayloadLimitx = [400, 400], [0, RangeArr[0]]
    PayRange, axis = plt.subplots()

    if RangeArr[-1] > 1000:
        maxRange_x = [1000, RangeArr[-1]]
        maxRange_y = [minpayload+availableSpace, minpayload+availableSpace]
        axis.plot(maxRange_x, maxRange_y, 'b-')

    axis.plot(RangeArr, payloadmass, 'b-')
    axis.plot(RangeLimitx, RangeLimity, label='Range requirement')
    axis.plot(PayloadLimitx, PayloadLimity, label='Payload requirement')
    plt.legend(loc='best')
    plt.ylim(390, RangeLimity[-1]+10)
    plt.xlim(left=750)
    plt.xlabel('Range [km]')
    plt.ylabel('Payload [kg]')
    plt.savefig('Tiltrotor-PayloadRangeDiagram.pdf')
    plt.show()


#Climb Performance
def AircraftClimbPerf(m_bat_cr, P_dens_cr, M_design, R, V_cr):
    W = M_design * const['gravityMars']
    P_climb_ac = m_bat_cr * P_dens_cr
    #Climb angle in cruise
    T_max = (P_climb_ac / (np.sqrt(1 / (2 * const['airDensity'] * np.pi * R**2))))**(2/3)
    gamma = np.arctan(T_max/W - const['cd']/const['cl'])
    ROC_cr = np.sin(gamma)*V_cr
    print(f'ROC in aircraft configuration: {ROC_cr} [m/s]')
    return ROC_cr


def RotorClimbPerf(MTOM, R, n_rotors):
    MTOW = MTOM * 2.205
    A = np.pi*R**2
    T_max = 1.1 * MTOM * const['gravityMars'] / 4
    T_hover = MTOM * const['gravityMars'] / 4
    V_c = 0*3.28084  # Climb Velocity
    v_1max = np.sqrt(T_max/(2*const['airDensity']*A))  # m/s
    v_1hover = np.sqrt(T_hover/(2*const['airDensity']*A))  # m/s

    DhpActual = n_rotors * (T_max*v_1max - T_hover*v_1hover) / 745.7
    Dhp = 0
    v_1hover *= 3.28084  # Convert to ft/s
    while Dhp <= DhpActual:
        V_c += 0.1
        Dhp = MTOW/550 * (V_c/2 + np.sqrt((V_c/2)**2 + v_1hover**2) - v_1hover)
    ROC_vert = V_c / 3.28084
    print(f'ROC in rotorcraft configuration: {ROC_vert} [m/s]')

    return ROC_vert
