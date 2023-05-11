import numpy as np
import scipy.integrate
from scipy.interpolate import InterpolatedUnivariateSpline

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

    #TEST
def RadiusMass(M, N_rotors):
    gm=3.721
    T_min = M*gm /N_rotors
    N_blades = 6
    b=N_blades
    T=0
    R=1
    while T<T_min:
        R+=0.1
        c=R/20
        theta_tip = 6*np.pi/180
        x0 = 0.1
        V_m = 220
        v_tip = 0.92*V_m
        omega = v_tip/R
        rho = 0.01
        n_elements = 10
        a0 = 6
        alpha0 = 2*np.pi/180
        A = np.pi*R**2
        r2R = np.arange(1, n_elements+1)/n_elements
        c2R = c/R
        M_local = (r2R)*(omega*R/V_m)
        a = a0/(np.sqrt(1-M_local**2))
        Dtheta = theta_tip/r2R
        theta0 = -min(Dtheta-np.ones(n_elements,)*alpha0)+theta_tip
        theta = theta0+Dtheta-alpha0
        v12Omegar = a*b*c2R/(16*np.pi*r2R)*(-1+np.sqrt(1+(32*np.pi*theta*r2R)/(a*b*c2R)))
        alpha = theta-np.arctan(v12Omegar)
        alpha = 6*np.pi/180 * np.ones(n_elements)
        cl = a*alpha
        cd = 0.025
        print(alpha*180/np.pi)
        DctDr2R = b*r2R**2*c2R*cl/(2*np.pi)

        funCT = InterpolatedUnivariateSpline(r2R, DctDr2R)
        Ct_lossless = funCT.integral(0, 1)
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

        DL = Ct*rho*(omega*R)**2

        Ct2sigma = Ct/(b*c/(np.pi*R))

        Cq = (CQ_profile+CQ_induced+DCQ_I)*.93

        T=rho*A*(omega*R)**2*Ct
        power = rho*A*(omega*R)**3*Cq

    return R, T, power#*0.00134

print(RadiusMass(3000, 4))

def twist(theta):
    alpha =  theta - np.arctan((theta))