# J109-1 verification: variable capacitor + two diodes charge pump
import sympy as sp
import numpy as np
C,E,d = sp.symbols('C E d', positive=True)
n = sp.symbols('n', integer=True, positive=True)
C0 = 2*C; k = 3
# symbolic recursion
Vn = sp.Function('V')
V1 = (C*E + C0*E)/(C/k + C0)
print('V1 =', sp.simplify(V1))
Vlim = sp.solve(sp.Eq(sp.Symbol('X')*(C/k+C0), C*E + C0*sp.Symbol('X')), sp.Symbol('X'))[0]
print('limit =', Vlim)
Vgen = 3*E - 2*E*sp.Rational(6,7)**n
print('V_n formula check n=1:', sp.simplify(Vgen.subs(n,1)-V1))
Vp = sp.Symbol('Vp')
print('recursion check:', sp.simplify(((C*E + C0*Vp)/(C/k+C0)).subs(Vp, Vgen) - Vgen.subs(n,n+1)))
# switching gaps
print('D2 starts at x =', 'd*V_{n-1}/E ; n=2 ->', sp.simplify(d*V1/E))
print('D1 starts (narrowing) at x =', 'k d E/V_n ; n=1 ->', sp.simplify(k*d*E/V1))
# work in widening of cycle 1 (energy method)
W1w = sp.Rational(1,2)*(C/k+C0)*V1**2 - sp.Rational(1,2)*(C+C0)*E**2
print('W widen (cycle1) =', sp.simplify(W1w))
# net work over cycle n
Vprev = Vgen.subs(n,n-1); Vcur = Vgen
Wn = sp.simplify(sp.Rational(1,2)*C0*(Vcur**2-Vprev**2) - E*C0*(Vcur-Vprev))
print('W_n =', sp.factor(Wn), ' W_1 =', sp.simplify(Wn.subs(n,1)))
Wtot = sp.Rational(1,2)*C0*((3*E)**2 - E**2) - E*C0*(3*E-E)
print('total W =', sp.simplify(Wtot), ' sum series =', sp.simplify(sp.summation(Wn,(n,1,sp.oo))))

# ---- quasi-static numerical simulation with ideal diodes ----
Cn, En, dn = 1.0, 1.0, 1.0
C0n = 2*Cn
def Cx(x): return Cn*dn/x
Qx = Cn*En; Q0 = C0n*En   # after initial charging: both at E
Wext = 0.0; Wbat = 0.0
Vs = []; sw = []
xs_up = np.linspace(dn, 3*dn, 40001); xs_dn = xs_up[::-1]
def settle(x, Qx, Q0, Wbat):
    c = Cx(x); VP = Qx/c; VQ = Q0/C0n
    if VP > VQ + 1e-15:          # D2 conducts
        V = (Qx+Q0)/(c+C0n); Qx = c*V; Q0 = C0n*V; VP = V
    if VP < En - 1e-15:          # D1 conducts
        dq = c*En - Qx; Wbat += En*dq; Qx = c*En
    return Qx, Q0, Wbat
for cyc in range(1, 31):
    d2on = None; d1on = None
    for path in (xs_up, xs_dn):
        for i in range(1, len(path)):
            x0, x1 = path[i-1], path[i]
            F = Qx**2/(2*Cn*dn)          # attractive force magnitude (charge-const)
            Wext += F*(x1-x0)
            Q0old = Q0; Qxold = Qx
            Qx, Q0, Wbat = settle(x1, Qx, Q0, Wbat)
            if path is xs_up and d2on is None and Q0 > Q0old + 1e-12: d2on = x0
            if path is xs_dn and d1on is None and Qx > Qxold + 1e-12: d1on = x0
    Vs.append(Q0/C0n); sw.append((d2on, d1on))
    if cyc == 1: W1 = Wext
print('numeric V1..V3:', Vs[:3], ' formula:', [float(Vgen.subs({n:i,E:1})) for i in (1,2,3)])
print('switch gaps cycle1, cycle2:', sw[0], sw[1], ' formula: D2 at 1, 9/7=%.4f ; D1 at 21/9=%.4f, %.4f' % (9/7, 3*7/9, 3/float(Vgen.subs({n:2,E:1}))))
print('numeric W cycle1 =', W1, ' formula =', float(Wn.subs({n:1,C:1,E:1})))
print('numeric W after 30 cycles =', Wext, ' formula total =', float(Wtot.subs({C:1,E:1})), '(approaches)')
U = 0.5*C0n*Vs[-1]**2 - 0.5*C0n*En**2
print('energy balance: Wext + Wbat =', Wext+Wbat, ' dU_res =', U)
