# J107-1 verification: charged disk shuttling between parallel plates
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
e0,a,V,d,m,g = sp.symbols('epsilon0 a V d m g', positive=True)
E = V/d
q = e0*a*E                       # charge on disk lying on B (surface charge density e0*E)
F_contact = sp.Rational(1,2)*e0*E**2*a   # force on surface charge = (1/2) sigma E
VL = sp.solve(sp.Eq(F_contact, m*g), V)[0]
print('q =', q, ' F_contact =', F_contact, ' V_L =', sp.simplify(VL))
# (3) V = sqrt2 VL
Vs = sp.sqrt(2)*VL
qs = sp.simplify(q.subs(V,Vs)); Fdet = sp.simplify((q*E).subs(V,Vs))
print('q at V=sqrt2 VL:', qs, ' detached force qE =', Fdet, ' /mg =', sp.simplify(Fdet/(m*g)))
a_up = sp.simplify((Fdet - m*g)/m); a_dn = sp.simplify((Fdet + m*g)/m)
print('a_up =', a_up, ' a_down =', a_dn)
T = sp.sqrt(2*d/a_up) + sp.sqrt(2*d/a_dn)
print('T =', sp.simplify(T), ' = sqrt(2d/g)*(1/sqrt3+1/sqrt5)?', sp.simplify(T - sp.sqrt(2*d/g)*(1/sp.sqrt(3)+1/sp.sqrt(5))))
I = 2*qs/T
print('I =', sp.simplify(I))
Iexp = 4*sp.sqrt(e0*a*m*g)/(sp.sqrt(2*d/g)*(1/sp.sqrt(3)+1/sp.sqrt(5)))
print('I check:', sp.simplify(I-Iexp))
# energy: KE at A, KE at B, battery work
qV = sp.simplify(qs*Vs); print('qV =', qV)
KEup = sp.simplify(m*a_up*d); KEdn = sp.simplify(m*a_dn*d)
print('KE at A =', KEup, ' KE at B =', KEdn, ' sum =', sp.simplify(KEup+KEdn), ' battery 2qV =', sp.simplify(2*qV))
# numerical check with event switching (state: on-B, flying up, on-A, flying down)
vals = {e0:8.85e-12, a:1e-4, m:2e-6, g:9.8, d:0.01}
VLn = float(VL.subs(vals)); Vn = np.sqrt(2)*VLn; dn = vals[d]
qn = float((q).subs(vals).subs(V,Vn)); En = Vn/dn; mn = vals[m]; gn = vals[g]
def run(sign):
    # sign=+1: disk +q flying from B (y=0) up to A (y=d); sign=-1: -q from A down
    f = lambda t,s: [s[1], (sign*qn*En)/mn - gn]
    ev = (lambda t,s: s[0]-dn) if sign>0 else (lambda t,s: s[0])
    ev.terminal=True
    y0 = [0.0,0.0] if sign>0 else [dn,0.0]
    sol = solve_ivp(f,[0,1],y0,events=ev,rtol=1e-11,atol=1e-14)
    return sol.t_events[0][0], sol.y_events[0][0][1]
t1,v1 = run(+1); t2,v2 = run(-1)
Tn = t1+t2; Tf = float(T.subs(vals))
print('numeric T =', Tn, ' formula =', Tf)
print('numeric KE sum =', 0.5*mn*(v1**2+v2**2), ' 2qV =', 2*qn*Vn)
print('contact force at B vs mg (must exceed):', 0.5*8.85e-12*En**2*1e-4, mn*gn)
