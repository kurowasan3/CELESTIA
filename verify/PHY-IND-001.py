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
# (4) general efficiency eta (u = -v here, left positive)
u1=-s[v1]; u2=-s[v2]
eta=sp.simplify(m2*u2/(m1*u1))
eta_claim=m2*(R*m1-(R+r)*m2)/(m1*((R+r)*m1-R*m2))
print('eta matches claim:',sp.simplify(eta-eta_claim)==0)
# eta < m2/m1 : numerator - denominator of (m1/m2)eta
print('num-den =',sp.factor((R*m1-(R+r)*m2)-((R+r)*m1-R*m2)))
# random check of inequality under lift condition
rng=np.random.default_rng(1)
ok=True
for _ in range(2000):
    Rv,rv,m2v=rng.uniform(0.1,5,3); m1v=m2v*(1+rv/Rv)*rng.uniform(1.001,5)
    e=float(eta_claim.subs({R:Rv,r:rv,m1:m1v,m2:m2v}))
    ok&= (0<e<m2v/m1v)
print('0<eta<m2/m1 for random lift-condition samples:',ok)
# r->0 limit
print('r->0: eta ->',sp.simplify(sp.limit(eta_claim,r,0)),' u1->',sp.factor(sp.limit(u1,r,0)),' u2->',sp.factor(sp.limit(u2,r,0)))
heatR=(R*s[U]**2/R**2)  # U^2/R
print('r->0: R heat / input =',sp.simplify(sp.limit((s[U]**2/R)/(m1*g*u1),r,0)))
print('u1-u2 =',sp.factor(u1-u2))
assert sp.simplify(eta-eta_claim)==0 and ok
