import numpy as np

Wcrew = 250 * 2.205 #kg
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