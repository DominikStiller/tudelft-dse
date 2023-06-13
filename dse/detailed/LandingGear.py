import numpy as np

maximumTirePressure = 2.8 #kg/cm^2
W = 2700
g = 3.71
S = 1
W2S = W*g / S
def TireSizing():
    pmax = 170000
    pTire =1

def LongitudinalPosition(xcg, ycg, zcg, b, ztip, lt):
    Pn = 0.08*W*g
    Pm = 0.92*W*g
    zlg = 1.5*zcg
    lm = lt - zlg/np.tan(np.radians(2.8))
    ln = (W*g / Pn - 1)*lm
    phi = np.radians(55)
    ymlgTIPOVER = (ln+lm)/(np.sqrt((ln**2 * np.tan(phi)**2)/(zcg**2) - 1))
    ymlgTIPCLEARANCE = b/2 - ztip/np.tan(np.degrees(5))

    print(ln, lm, zlg, max(ymlgTIPCLEARANCE, ymlgTIPOVER))

LongitudinalPosition(0, 0, 0.3, 42.69, 3, 10)


