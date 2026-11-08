# T20082-1 battery E - resistor R - capacitors A,B in series; neon lamp across A; flashes instantaneous
import numpy as np
from fractions import Fraction as Fr

C = 1.0          # results in units of C (energies in C*V^2)
CA, CB = 1*C, 2*C
Von, Voff, E = 90.0, 60.0, 280.0
Rres = 1.0
dt = 1e-4

# state: QA (charge on top plate of A), QB (charge on top plate of B)
# node (A bottom + B top) charge q = QB - QA ; top terminal potential relative to bottom of B: Vt = QA/CA + QB/CB
QA = QB = 0.0
heat = 0.0; work = 0.0; lamp = 0.0; flashes = 0
t = 0.0
first_flash_VB = None
while t < 60.0:   # RC_series = R*CA*CB/(CA+CB) = 2/3 ; 60 >> that
    Vt = QA/CA + QB/CB
    i = (E - Vt)/Rres
    # RK2 for charge through R: same current charges QA and QB (series)
    Vt_mid = (QA + i*dt/2)/CA + (QB + i*dt/2)/CB
    i_mid = (E - Vt_mid)/Rres
    QA += i_mid*dt; QB += i_mid*dt
    heat += i_mid**2*Rres*dt; work += E*i_mid*dt
    t += dt
    if QA/CA >= Von:
        if first_flash_VB is None:
            first_flash_VB = QB/CB
        # instantaneous flash: lamp moves charge from top plate of A to node until V_A = Voff ; B untouched
        U_before = 0.5*QA**2/CA
        QA = CA*Voff
        lamp += U_before - 0.5*QA**2/CA
        flashes += 1
VA = QA/CA; VB = QB/CB
stored = 0.5*CA*VA**2 + 0.5*CB*VB**2
print('flashes', flashes, ' V_B at first flash', first_flash_VB)
print('final VA, VB', VA, VB, ' QB final', QB)
print('battery work', work, ' E*QB', E*QB)
print('lamp total', lamp, ' per flash', lamp/flashes)
print('stored', stored, ' heat sim', heat, ' heat from balance', work - stored - lamp)

# exact (fractions)
CA, CB, Von, Voff, E = Fr(1), Fr(2), Fr(90), Fr(60), Fr(280)
dq = CA*(Von - Voff)
N = 0
while (CB*E - N*dq)/(CA + CB) > Von:
    N += 1
VAf = (CB*E - N*dq)/(CA + CB); VBf = E - VAf
W = E*CB*VBf; St = CA*VAf**2/2 + CB*VBf**2/2; Lp = N*CA*(Von**2 - Voff**2)/2
print('exact N', N, ' VAf', VAf, float(VAf), ' VBf', VBf, ' W', W, float(W), ' stored', St, float(St), ' lamp', Lp, ' heat', W - St - Lp, float(W - St - Lp))
