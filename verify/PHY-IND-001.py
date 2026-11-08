# N36-1: two rods P1,P2 (resistance r each, mass M) on horizontal rails (spacing l, no resistance), resistor R between rails.
# B vertical. P1 pulled toward -x by m1 (string over pulley at left), P2 pulled toward +x by m2 (pulley at right).
# conventions: rod current i_k (+y : X->Y) -> force i_k*l*B (+x).  emf of rod moving v (+x): V_Y-V_X = -vBl - r i
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
B,l,r,R,g,m1,m2,M=sp.symbols('B l r R g m1 m2 M',positive=True)
v1,v2,i1,i2,U=sp.symbols('v1 v2 i1 i2 U')
eqs=[sp.Eq(U,-v1*B*l-r*i1),sp.Eq(U,-v2*B*l-r*i2),sp.Eq(i1+i2,U/R),
     sp.Eq(i1*l*B,m1*g),sp.Eq(i2*l*B,-m2*g)]   # steady force balance
s=sp.solve(eqs,[v1,v2,i1,i2,U],dict=True)[0]
for k in s: print(k,'=',sp.factor(s[k]))
# lift condition: v2<0
print('v2<0 <=>', sp.factor(sp.numer(sp.together(-s[v2]))))
# energy balance: m1 g |v1| = m2 g |v2| + Joule
J=r*s[i1]**2+r*s[i2]**2+s[U]**2/R
print('energy check:',sp.simplify(m1*g*(-s[v1]) - (m2*g*(-s[v2])) - J))
print('Joule=',sp.factor(J))
# numeric ODE check
vals={B:0.5,l:0.4,r:0.2,R:0.3,g:9.8,m1:0.03,m2:0.01,M:0.02}
Bn,ln,rn,Rn,gn,m1n,m2n,Mn=[vals[x] for x in (B,l,r,R,g,m1,m2,M)]
def rhs(t,y):
    V1,V2=y
    # solve circuit given velocities: unknowns i1,i2,U
    A=np.array([[rn,0,1],[0,rn,1],[1,1,-1/Rn]]); b=np.array([-V1*Bn*ln,-V2*Bn*ln,0])
    I1,I2,Uu=np.linalg.solve(A,b)
    a1=(I1*ln*Bn - m1n*gn)/(Mn+m1n); a2=(I2*ln*Bn + m2n*gn)/(Mn+m2n)
    return [a1,a2]
so=solve_ivp(rhs,[0,20],[0,0],rtol=1e-10,atol=1e-12)
print('ODE final v1,v2=',so.y[:,-1],' formula:',float(s[v1].subs(vals)),float(s[v2].subs(vals)))
# case m1=3m2, r=R : efficiency
sub={m1:3*m2,r:R}
eff=sp.simplify((m2*g*(-s[v2]))/(m1*g*(-s[v1])).subs(sub)) if False else sp.simplify(((m2*(-s[v2]))/(m1*(-s[v1]))).subs(sub))
print('v1,v2 at m1=3m2,r=R:',sp.simplify(s[v1].subs(sub)),sp.simplify(s[v2].subs(sub)),' efficiency=',eff)
