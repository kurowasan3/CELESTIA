# K20011-1 verification: dynamo with coil inductance -> current-source behaviour
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
N,Ph,L,r,w,r1,r2,R,a,b=sp.symbols('N Phi0 L r omega r1 r2 R a b',positive=True)
p=4*w; E0=2*N*Ph*p; X=2*L*p
I0=E0/sp.sqrt(r**2+X**2); P=sp.Rational(1,2)*I0**2*r
print('I0=',sp.simplify(I0),' P=',sp.simplify(P))
print('lim I0 w->oo:',sp.limit(I0,w,sp.oo),' lim P:',sp.limit(P,w,sp.oo))
print('half sat w:',sp.solve(sp.Eq(P,sp.limit(P,w,sp.oo)/2),w))
Pr=lambda rr: P.subs(r,rr)
print('crossover w:',sp.solve(sp.Eq(Pr(r1),Pr(r2)),w))
print('v for w*: ',sp.simplify(R*a/b*sp.sqrt(r1*r2)/(8*L)))
# ODE check: 2L di/dt + r i = emf(t) where emf = -d/dt(2N Phi0 cos pt)
vals=dict(N=100,Ph=1e-5,L=0.02,r=10.0,w=150.0)
pp=4*vals['w']
def rhs(t,y): return [(2*vals['N']*vals['Ph']*pp*np.sin(pp*t)-vals['r']*y[0])/(2*vals['L'])]
sol=solve_ivp(rhs,[0,0.3],[0],max_step=1e-5,rtol=1e-9)
amp=np.max(np.abs(sol.y[0][sol.t>0.2]))
print('num amp',amp,' formula',float(I0.subs({N:100,Ph:1e-5,L:0.02,r:10.0,w:150.0})),' sat',100*1e-5/0.02)
# flux linkage check at high omega: 2N Phi + 2L i ~ const
vals['w']=5000.0; pp=4*vals['w']
sol=solve_ivp(rhs,[0,0.06],[0],max_step=2e-7,rtol=1e-9)
t=sol.t; i=sol.y[0]; m=t>0.05
lam=2*vals['N']*vals['Ph']*np.cos(pp*t)+2*vals['L']*i
print('high-w flux linkage variation (rel to 2NPhi0):',(lam[m].max()-lam[m].min())/(2*vals['N']*vals['Ph']))
