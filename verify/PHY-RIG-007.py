# T20021-1 verification: pipe with hidden balls m, 2m as a compound pendulum
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

m, g, L, l1, l2, th, x = sp.symbols('m g L l1 l2 theta x', positive=True)
M = 3*m
lc = (m*l1 + 2*m*l2)/M
I_A = m*l1**2 + 2*m*l2**2
# small oscillation via energy: E = 1/2 I thd^2 + M g lc (1-cos th)
Leq = sp.simplify(I_A/(M*lc))
print('equivalent length about A =', Leq)
# inverse problem: given l (=lc) and lam (=Leq)
l, lam = sp.symbols('l lambda', positive=True)
sols = sp.solve([sp.Eq(l1 + 2*l2, 3*l), sp.Eq(l1**2 + 2*l2**2, 3*l*lam)], [l1, l2], dict=True)
for s in sols:
    print('solution:', sp.simplify(s[l1]), ',', sp.simplify(s[l2]))
# moment about CM, radius of gyration
k2 = sp.simplify((I_A - M*lc**2)/M)
print('k^2 =', sp.factor(k2), '  == l(lam-l)?', sp.simplify(k2.subs({l1: sols[0][l1], l2: sols[0][l2]}) - l*(lam - l)) == 0)
# period about pivot at distance s from CM: Leq(s) = s + k^2/s
s = sp.symbols('s', positive=True)
K2 = sp.symbols('K2', positive=True)
Leq_s = s + K2/s
# other pivot with same Leq on the other side: s' = K2/s ; distance between = s + K2/s = Leq
print('s_conj = K2/s ; distance between pivots =', sp.simplify(s + K2/s), ' = Leq(s)')
smin = sp.solve(sp.diff(Leq_s, s), s)
print('min at s =', smin, ' Leq_min =', sp.simplify(Leq_s.subs(s, smin[0])))
# pivot at B: equivalent length
lamB = sp.simplify(((L - l)**2 + l*(lam - l))/(L - l))
print('lambda_B =', lamB)
# numeric check: l1=0.2, l2=0.7, L=1 -> simulate pendulum period about A and about conjugate point
gv = 9.8; l1v, l2v, Lv = 0.2, 0.7, 1.0
lcv = (l1v + 2*l2v)/3
lamv = (l1v**2 + 2*l2v**2)/(3*lcv)
def period(pivot):
    d1, d2 = l1v - pivot, l2v - pivot      # signed positions from pivot along pipe
    I = d1**2 + 2*d2**2
    mom = (d1 + 2*d2)                       # lever arm sum (sign: CM below pivot if >0)
    f = lambda t, y: [y[1], -gv*abs(mom)/I*np.sin(y[0])]
    sol = solve_ivp(f, [0, 10], [0.01, 0], events=lambda t, y: y[1], rtol=1e-11, atol=1e-13, dense_output=True)
    te = [t for t in sol.t_events[0] if t > 1e-6]
    return 2*te[0]
TA = period(0.0)
print('T_A numeric', TA, ' theory', 2*np.pi*np.sqrt(lamv/gv))
P = lamv   # conjugate pivot at distance lam from A (beyond CM)
print('pivot P at', P, '(within pipe:', P <= Lv, ')  T_P numeric', period(P))
print('recovered positions:', [(float(sp.N(sv[l1].subs({l: lcv, lam: lamv}))), float(sp.N(sv[l2].subs({l: lcv, lam: lamv})))) for sv in sols])
# both solutions give identical periods about any pivot
alt = [float(sp.N(sv[l1].subs({l: lcv, lam: lamv}))) for sv in sols]
alt2 = [float(sp.N(sv[l2].subs({l: lcv, lam: lamv}))) for sv in sols]
for p in [0.0, 0.5, 1.0]:
    out = []
    for A1, A2 in zip(alt, alt2):
        d1, d2 = A1 - p, A2 - p
        out.append((d1**2 + 2*d2**2)/abs(d1 + 2*d2))
    print('pivot', p, 'equiv lengths of both solutions', out)
