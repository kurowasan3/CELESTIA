# T20042-1 ion in crossed E (+z) and B (+x): turning points via two first integrals
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
q,m,E,B,V,h,z,k=sp.symbols('q m E B V h z kappa',positive=True)
v0=sp.sqrt(2*q*V/m)
# first integrals: vy = v0 + qBz/m ; (1/2)m(vy^2+vz^2) = qV + qEz
vy=v0+q*B*z/m
eq=sp.Eq(sp.Rational(1,2)*m*vy**2, q*V+q*E*z)   # vz=0 turning point
zs=sp.solve(eq,z); print('turning z:',zs)
zstar=[s for s in zs if s!=0][0]
target=2*m*(E-v0*B)/(q*B**2)
print('z* matches:',sp.simplify(zstar-target)==0)
# straight going q/m
ks=sp.solve(sp.Eq(sp.sqrt(2*k*V)*B,E),k); print('straight kappa:',ks)
# upper plate condition z* >= h  with kappa=q/m
zk=2*(E-B*sp.sqrt(2*k*V))/(k*B**2)
s=sp.symbols('s',positive=True)
roots=sp.solve(sp.Eq(zk.subs(k,s**2),h),s); print('roots s:',roots)
kup=2*(sp.sqrt(V+h*E)-sp.sqrt(V))**2/(h**2*B**2)
print('kappa_upper matches:',[sp.simplify(r**2-kup) for r in roots])
# deepest descent over kappa
depth=-zk
sc=sp.solve(sp.diff(depth.subs(k,s**2),s),s); print('critical s:',sc)
print('max depth:',sp.simplify(depth.subs(k,sc[0]**2)), ' kappa*/kappa_straight =',sp.simplify(sc[0]**2/ks[0]))
# numeric ODE check
qv,mv,Ev,Bv,Vv=1.0,1.0,2.0,1.0,1.0
def run(kap,T=60):
    v0n=np.sqrt(2*kap*Vv)
    f=lambda t,Y:[0,Y[4],Y[5],0,kap*Bv*Y[5],kap*(Ev-Bv*Y[4])]  # x,y,z,vx,vy,vz ; F=q(E+v x B)
    sol=solve_ivp(f,[0,T],[0,0,0,0,v0n,0],max_step=0.002,rtol=1e-10,atol=1e-12)
    return sol
for kap in [0.5,2.0,8.0,20.0]:
    sol=run(kap); zmin,zmax=sol.y[2].min(),sol.y[2].max()
    pred=float(zk.subs({E:Ev,B:Bv,V:Vv,k:kap}))
    print(kap,'zmin,zmax',zmin,zmax,'pred',pred)
    assert abs((zmax if pred>0 else zmin)-pred)<1e-4
# max depth numeric
res=minimize_scalar(lambda kap: run(kap,30).y[2].min(),bounds=(0.5,40),method='bounded')
print('numeric deepest kappa',res.x,'depth',-res.fun,' predicted kappa',2*Ev**2/(Vv*Bv**2),'depth',Vv/Ev)
assert abs(-res.fun-Vv/Ev)<1e-3
# upper plate threshold numeric
hv=0.7
kth=float(kup.subs({E:Ev,B:Bv,V:Vv,h:hv}))
print('k threshold',kth,' zmax at kth',run(kth).y[2].max())
print('OK')
