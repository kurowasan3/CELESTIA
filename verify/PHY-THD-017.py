# T20023-1 verification: gas pushing liquid into an evacuated upper chamber (P proportional to V)
import sympy as sp
rho, g, S, l, T0, R, Cv, x, H, n = sp.symbols('rho g S ell T0 R C_v x H n', positive=True)
P1 = rho*g*l                    # initial: liquid in tube stands l above B surface, A vacuum
V1 = S*l
nR = P1*V1/T0
# phase 1 (isochoric) ends when P = rho g * 2l (tube top = A floor at 2l above B surface)
Pa = rho*g*2*l
Ta = sp.simplify(Pa*V1/nR)
print('T at switch =', Ta)
# phase 2: B level drops x, A level rises x (same cross-section S); head = 2l + 2x ; V = S(l + x)
V = S*(l + x)
P = rho*g*(2*l + 2*x)
print('P/V =', sp.simplify(P/V), ' (const => P propto V)')
T = sp.simplify(P*V/nR)
print('T(V) =', sp.factor(T))
# molar heat along process: dQ = n Cv dT + P dV
dQdx = n*Cv*sp.diff(T, x) + P*sp.diff(V, x)
C = sp.simplify((dQdx/(n*sp.diff(T, x))).subs(n, nR/R))
print('molar heat in phase 2 =', C)
r = sp.symbols('r')
Cv_sol = sp.solve(sp.Eq(Cv/(Cv + R/2), sp.Rational(5, 6)), Cv)
print('Cv from ratio 5/6 =', Cv_sol)
# total heat until V = 2 S l (x = l)
Q1 = (nR/R)*Cv*(Ta - T0)
Q2 = sp.integrate(dQdx.subs(n, nR/R), (x, 0, l))
Qtot = sp.simplify((Q1 + Q2).subs(Cv, Cv_sol[0]))
print('Q1+Q2 =', sp.factor(Qtot))
# work = gain of PE of liquid (A is vacuum)
W = sp.integrate(P*S, (x, 0, l))
# PE: slab of thickness x moved from B surface region (heights -x..0 rel. B initial surface) to A (2l..2l+x)
y = sp.symbols('y')
PE = rho*g*S*(sp.integrate(y, (y, 2*l, 3*l)) - sp.integrate(y, (y, -l, 0)))
print('W =', sp.simplify(W), '  dPE =', sp.simplify(PE), ' equal:', sp.simplify(W - PE) == 0)
# general H
Pg = rho*g*(H + 2*x); Vg = S*(l + x)
nRg = P1*V1/T0
Tg = Pg*Vg/nRg
Cg = sp.simplify((Cv + R*Pg*sp.diff(Vg, x)/(sp.diff(Pg*Vg, x))).subs(x, 0))
print('general H: molar heat just after switch =', sp.simplify(Cg), ' ratio r2/r1 =', sp.simplify(Cv/Cg))
print('  H=2l check:', sp.simplify(Cg.subs(H, 2*l) - (Cv + R/2)) == 0)
# check temperatures positive & plausibility: Tend = 8 T0
print('T at x=l:', sp.simplify(T.subs(x, l)))
# numeric check via small-step integration of first law with constant heater power
import numpy as np
rv, gv, Sv, lv, T0v, Rv = 1.36e4, 9.8, 1e-3, 0.1, 300.0, 8.31
Cvv = 2.5*Rv
P1v = rv*gv*lv; nv = P1v*Sv*lv/(Rv*T0v)
Ph = 1.0  # W
# phase 1 rate
r1 = Ph/(nv*Cvv)
# phase 2: simulate dx/dt from dQ = n Cv dT + P dV
xs = 0.0; Tn = 2*T0v; dt = 1e-3; Ts = []
for k in range(2000):
    Pn = rv*gv*(2*lv + 2*xs)
    dTdx = (rv*gv*Sv*(2*lv + 2*xs) + rv*gv*2*Sv*(lv + xs))/(nv*Rv)
    dQdx_ = nv*Cvv*dTdx + Pn*Sv
    dx = Ph*dt/dQdx_
    xs += dx; Ts.append(Tn); Tn = rv*gv*(2*lv + 2*xs)*Sv*(lv + xs)/(nv*Rv)
r2 = (Ts[100] - Ts[0])/(100*dt)
print('numeric r2/r1 =', r2/r1, ' expected 5/6 =', 5/6)
