# J117-1 verification: two LEDs (piecewise-linear) in parallel with common series resistor
import sympy as sp
import numpy as np
from scipy.optimize import brentq
V1t, V2t, r, R = sp.Rational(18,10), sp.Rational(24,10), 10, 20
E, V, Rs = sp.symbols('E V R_s', positive=True)
# (1) LED2 threshold
I1 = (V2t-V1t)/r; E2 = V2t + R*I1
print('(1) E at LED2 on =', E2, float(E2), ' I1 then =', I1)
# (2) E=6.0
sol = sp.solve(sp.Eq((6-V)/R, (V-V1t)/r + (V-V2t)/r), V)[0]
i1 = (sol-V1t)/r; i2 = (sol-V2t)/r
print('(2) VP =', sol, float(sol), ' I1 =', float(i1), ' I2 =', float(i2), ' ratio', float(i2/i1))
# (3) ratio bound with I1max = 0.20 A
VP = V1t + r*sp.Rational(2,10); print('(3) VP at I1=0.2:', VP, ' I2=', (VP-V2t)/r, ' ratio', (VP-V2t)/r/sp.Rational(2,10))
EmaxP = VP + R*(sp.Rational(2,10)+(VP-V2t)/r); print('   E at that point =', EmaxP)
# general: ratio = (V-2.4)/(V-1.8) < 1 always
# (4) series Rs in LED1 branch, equal currents at E=6.0
I = sp.symbols('I', positive=True)
s = sp.solve([sp.Eq(V2t + r*I, V1t + (r+Rs)*I), sp.Eq(6 - (V2t + r*I), R*2*I)], [I, Rs], dict=True)[0]
print('(4) I =', s[I], float(s[I]), ' Rs =', s[Rs], float(s[Rs]))
Rsn = float(s[Rs])
# (5) numeric: for E = 5 and 7, compare currents (piecewise model solved numerically)
def led(Vv, Vt, rr): return max(0.0, (Vv-Vt)/rr)
def currents(En, Rsv):
    # LED1 branch: series Rs+r, threshold 1.8 ; LED2 branch r, threshold 2.4
    f = lambda Vp: (En-Vp)/R - led(Vp,1.8,r+Rsv) - led(Vp,2.4,r)
    Vp = brentq(f, 0, En); return led(Vp,1.8,r+Rsv), led(Vp,2.4,r), Vp
for En in (3.6, 5.0, 6.0, 7.0):
    print('  E=%.1f, Rs=0: I1,I2,VP=' % En, currents(En,0.0), '  Rs=25/3:', currents(En,Rsn))
# E at which LED2 turns on with Rs inserted
E2s = 2.4 + R*(0.6/(r+Rsn)); print('(5) with Rs, LED2 on at E =', E2s)
