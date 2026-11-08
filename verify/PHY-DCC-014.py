# J120-1 verification: solar cell + C + R in parallel; light on then off
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
R,C,V0,t = sp.symbols('R C V_0 t', positive=True)
r = R; sP = 2*V0/R          # given: r = R, sP*R = 2 V0
# phase 1: current source sP into R||C
V1t = sP*R*(1-sp.exp(-t/(R*C)))
t1 = sp.solve(sp.Eq(V1t, V0), t)[0]; print('t1 =', sp.simplify(t1))
# phase 2: I = sP - (V-V0)/r
Vinf = sp.solve(sp.Eq(sP - (sp.Symbol('X')-V0)/r, sp.Symbol('X')/R), sp.Symbol('X'))[0]
print('V_inf =', Vinf, ' tau2 = C*rR/(r+R) =', sp.simplify(C*r*R/(r+R)))
Iinf = sp.simplify(sP - (Vinf-V0)/r); print('I_cell inf =', Iinf)
# light off: I = -(V-V0)/r for V>V0 ; else 0
Itjump = sp.simplify(-(Vinf-V0)/r); print('cell current just after off =', Itjump)
Vt = V0*R/(r+R); tau_a = C*r*R/(r+R)
ta = sp.solve(sp.Eq(Vt + (Vinf-Vt)*sp.exp(-t/tau_a), V0), t)[0]; print('time to reach V0 after off =', sp.simplify(ta))
Qtot = C*(Vinf-V0)
# Q_R - Q_cell = int (V/R - (V-V0)/r) dt = V0/R * ta (since r=R)
diff = V0/R*ta
QR = sp.simplify((Qtot+diff)/2); Qc = sp.simplify((Qtot-diff)/2)
print('Q_R =', QR, ' Q_cell =', Qc)
# direct integration check
Va = Vt + (Vinf-Vt)*sp.exp(-t/tau_a)
print('integral check Q_R:', sp.simplify(sp.integrate(Va/R,(t,0,ta)) - QR), ' Q_cell:', sp.simplify(sp.integrate((Va-V0)/r,(t,0,ta)) - Qc))

# numerical: piecewise ODE
Rn, Cn, V0n = 1.0, 1.0, 1.0; rn = Rn; sPn = 2*V0n/Rn
def cell(V, light):
    s = sPn if light else 0.0
    return s if V <= V0n else s - (V-V0n)/rn
toff = 20.0
def f(tt, y, light): return [(cell(y[0],light) - y[0]/Rn)/Cn, cell(y[0],light) if not light else 0.0, (y[0]/Rn) if not light else 0.0]
ev = lambda tt,y,l: y[0]-V0n
s1 = solve_ivp(lambda tt,y: f(tt,y,True), [0,toff], [0,0,0], events=lambda tt,y: y[0]-V0n, rtol=1e-11, atol=1e-13, max_step=1e-3)
print('numeric t1 =', s1.t_events[0][0], ' formula', float(t1.subs({R:1,C:1})))
print('numeric V before off =', s1.y[0,-1], ' formula', float(Vinf.subs(V0,1)))
evb = lambda tt,y: y[0]-V0n; evb.terminal=True; evb.direction=-1
s2 = solve_ivp(lambda tt,y: f(tt,y,False), [0,10], [s1.y[0,-1],0,0], events=evb, rtol=1e-11, atol=1e-13, max_step=1e-3)
print('numeric time to V0 after off =', s2.t_events[0][0], ' formula', float(ta.subs({R:1,C:1})))
print('numeric Q_R, Q_cell(absorbed) =', s2.y_events[0][0][2], -s2.y_events[0][0][1], ' formula', float(QR.subs({C:1,V0:1})), float(Qc.subs({C:1,V0:1})))
