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

# ---- (5) rewritten: energy balance of each charging interval between flashes ----
import sympy as sp
CAs,CBs,Es,Von_s,Voff_s,k=sp.symbols('C_A C_B E V_on V_off k',positive=True)
dV=Von_s-Voff_s; q=k*CAs*dV; QA=sp.symbols('Q_A')
U=QA**2/(2*CAs)+(q+QA)**2/(2*CBs)
dU=sp.simplify(U.subs(QA,CAs*Von_s)-U.subs(QA,CAs*Voff_s))
W=Es*CAs*dV
H=sp.simplify(W-dU)
H_ans=CAs*dV*(Es-(CAs+CBs)*(Von_s+Voff_s)/(2*CBs)-k*CAs*dV/CBs)
print('(5a) H_k - answer =',sp.simplify(H-H_ans))
assert sp.simplify(H-H_ans)==0
dec=sp.simplify(H_ans-H_ans.subs(k,k+1)); print('(5b) decrement',sp.factor(dec))
assert sp.simplify(dec-CAs**2*dV**2/CBs)==0
# numeric simulation (C=1, R=1): heat and battery work in each interval between flashes
CA_, CB_, Von_, Voff_, E_ = 1.0, 2.0, 90.0, 60.0, 280.0
QA_=QB_=0.0; h=w=0.0; t=0.0; dt=2e-5; marks=[]
while t<40.0:
    def cur(a,b): return (E_-(a/CA_+b/CB_))
    k1=cur(QA_,QB_); k2=cur(QA_+k1*dt/2,QB_+k1*dt/2); k3=cur(QA_+k2*dt/2,QB_+k2*dt/2); k4=cur(QA_+k3*dt,QB_+k3*dt)
    i=(k1+2*k2+2*k3+k4)/6
    QA_+=i*dt; QB_+=i*dt; h+=((k1**2+4*k2**2+k4**2)/6)*dt; w+=E_*i*dt; t+=dt
    if QA_/CA_>=Von_:
        marks.append((h,w)); QA_=CA_*Voff_
print('flashes',len(marks))
for kk in range(1,len(marks)):
    Hn=marks[kk][0]-marks[kk-1][0]; Wn=marks[kk][1]-marks[kk-1][1]
    Hf=float(H_ans.subs({CAs:CA_,CBs:CB_,Es:E_,Von_s:Von_,Voff_s:Voff_,k:kk}))
    print(f'  interval k={kk}: H num {Hn:.3f}  formula {Hf:.3f}   W num {Wn:.3f}  formula {E_*CA_*(Von_-Voff_):.1f}')
    assert abs(Hn-Hf)<2e-3*Hf+1 and abs(Wn-E_*CA_*(Von_-Voff_))<1e-3*Wn
print('CAP-013 (5) OK')
