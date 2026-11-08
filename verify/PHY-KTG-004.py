# T20063-1 Knudsen effusion (1/6 model), force on partition independent of its temperature
import sympy as sp
import random

R, Mm, NA, T1, T2, Tp, V1, V2, n, S, s = sp.symbols('R M N_A T_1 T_2 T_p V_1 V_2 n S s', positive=True)
mu = Mm/NA
v = lambda T: sp.sqrt(3*R*T/Mm)          # 1/6 model: all molecules speed v, (1/2)mu v^2 = (3/2)kT
N1, N2 = sp.symbols('N1 N2', positive=True)  # molecule numbers
nd1, nd2 = N1/V1, N2/V2
P1 = sp.Rational(1, 3)*nd1*mu*v(T1)**2
P2 = sp.Rational(1, 3)*nd2*mu*v(T2)**2
# check ideal-gas law consistency P V = (N/NA) R T
print('PV=nRT check', sp.simplify(P1*V1 - N1/NA*R*T1))
# steady state: equal number flux through hole
flux1 = sp.Rational(1, 6)*nd1*v(T1)*s
flux2 = sp.Rational(1, 6)*nd2*v(T2)*s
solN = sp.solve([sp.Eq(flux1, flux2), sp.Eq(N1 + N2, n*NA)], [N1, N2], dict=True)[0]
P1s = sp.simplify(P1.subs(solN)); P2s = sp.simplify(P2.subs(solN))
print('P1/P2', sp.simplify(P1s/P2s), ' vs sqrt(T1/T2)', sp.simplify(P1s/P2s - sp.sqrt(T1/T2)))
P1claim = n*R/(V1/T1 + V2/sp.sqrt(T1*T2))
print('P1 check', sp.simplify(P1s - P1claim))
# force on partition (area S, both faces), faces at Tp, molecules stick then re-emitted normally at v(Tp)
def face_pressure(nd, T):
    inc = sp.Rational(1, 6)*nd*v(T)       # number per area per time
    return inc*mu*v(T) + inc*mu*v(Tp)
F = sp.simplify((face_pressure(nd1, T1) - face_pressure(nd2, T2)).subs(solN)*S)
print('F - (P1-P2)S/2 =', sp.simplify(F - (P1s - P2s)*S/2))
# energy flux through hole (net) from 1 to 2
Edot = sp.simplify((flux1*sp.Rational(1, 2)*mu*v(T1)**2 - flux2*sp.Rational(1, 2)*mu*v(T2)**2).subs(solN))
Eclaim = sp.Rational(1, 4)*P1s*v(T1)*s*(T1 - T2)/T1
print('energy flux check', sp.simplify(Edot - Eclaim))
# random numeric spot check
for _ in range(3):
    vals = {R: 8.31, Mm: random.uniform(0.002, 0.05), NA: 6.02e23, T1: random.uniform(300, 600), T2: random.uniform(100, 300),
            Tp: random.uniform(50, 900), V1: random.uniform(1e-3, 1e-2), V2: random.uniform(1e-3, 1e-2), n: 1e-4, S: 1e-2, s: 1e-8}
    print(float(F.subs(vals)), float(((P1s - P2s)*S/2).subs(vals)))
