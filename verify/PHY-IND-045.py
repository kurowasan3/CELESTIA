# S0522-1 verification: massless rods + spring in B behave as inductor L=B^2 l^2/k
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp, quad
B,l,k,R,r,E,t,i=sp.symbols('B l k R r E t i',positive=True)
L=B**2*l**2/k
# (1) S2 open: E - B l x' = R i,  B l i = k x  (massless: force balance)
x=sp.Function('x')
ode=sp.Eq(B*l*(E-B*l*x(t).diff(t))/R, k*x(t))
xs=sp.dsolve(ode,ics={x(0):0}).rhs; print('x(t)=',sp.simplify(xs))
it=sp.simplify(k*xs/(B*l)); print('i(t)=',it,' check E=Ri+L di/dt:',sp.simplify(E-R*it-L*sp.diff(it,t)))
print("x'(0+)=",sp.simplify(sp.diff(xs,t).subs(t,0)))
# (2) S2 closed: N current B l x'/r gives force -B^2l^2x'/r
ode2=sp.Eq(B*l*(E-B*l*x(t).diff(t))/R - B**2*l**2*x(t).diff(t)/r, k*x(t))
xs2=sp.dsolve(ode2,ics={x(0):0}).rhs; xs2=sp.simplify(xs2); print('x2(t)=',xs2)
v0=sp.simplify(sp.diff(xs2,t).subs(t,0)); print("x2'(0+)=",v0)
i0=sp.simplify((E-B*l*v0)/R); ir0=sp.simplify(B*l*v0/r); print('i(0+)=',i0,' i_r(0+)=',ir0)
tau=sp.simplify(L*(R+r)/(R*r)); print('tau check:',sp.simplify(xs2-(B*l*E/(k*R))*(1-sp.exp(-t/tau))))
xinf=B*l*E/(k*R)
# (3) S1 to b: B l(-B l x')/R - B^2 l^2 x'/r = k x  -> heats
vals={B:0.5,l:0.4,k:3.0,R:2.0,r:5.0,E:6.0}
Bn,ln,kn,Rn,rn,En=[vals[s] for s in (B,l,k,R,r,E)]
xi=Bn*ln*En/(kn*Rn)
f=lambda tt,s: [-kn*s[0]/(Bn**2*ln**2*(1/Rn+1/rn))]
sol=solve_ivp(f,[0,60],[xi],rtol=1e-12,atol=1e-14,dense_output=True)
ts=np.linspace(0,60,600001); xx=sol.sol(ts)[0]; vv=-kn*xx/(Bn**2*ln**2*(1/Rn+1/rn))
Hr=np.trapezoid((Bn*ln*vv)**2/rn,ts); HR=np.trapezoid((Bn*ln*vv)**2/Rn,ts)
U=0.5*kn*xi**2
print('H_r num',Hr,' pred',U*Rn/(Rn+rn)); print('H_R num',HR,' pred',U*rn/(Rn+rn))
print('U =',sp.simplify(sp.Rational(1,2)*k*xinf**2), ' = 1/2 L (E/R)^2 ?',sp.simplify(sp.Rational(1,2)*k*xinf**2-L*(E/R)**2/2))
