# Q064-1: caterpillar on uniformly stretched string: (a) winch reels in at fixed position, (b) end moves away
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
t,L,v0,c,x0=sp.symbols('t L v0 c x0',positive=True)
X=sp.Function('X')
# case (a): dx/dt = v0 x/L - c
sol=sp.dsolve(sp.Eq(X(t).diff(t),v0*X(t)/L-c),X(t),ics={X(0):x0}); print(sol)
ta=sp.solve(sp.Eq(sol.rhs.subs(x0,c*L/(2*v0)),0),t); print('time (a), x0=cL/2v0:',ta)
# case (b): L(t)=L+v0 t, dx/dt = v0 x/(L+v0 t) - c
solb=sp.dsolve(sp.Eq(X(t).diff(t),v0*X(t)/(L+v0*t)-c),X(t),ics={X(0):x0}); print(sp.simplify(solb.rhs))
tb=sp.solve(sp.Eq(sp.simplify(solb.rhs),0),t); print('time (b):',tb)
# numbers: L=1 m, v0=0.01, c=0.001, x0=0.5
Ln,vn,cn,xn=1.0,0.01,0.001,0.5
print('(a) threshold x0 <',cn*Ln/vn)
ev=lambda tt,y:y[0]; ev.terminal=True
r=solve_ivp(lambda tt,y:[vn*y[0]/(Ln+vn*tt)-cn],[0,1e6],[xn],events=ev,rtol=1e-10,atol=1e-12)
print('(b) numeric time',r.t_events[0], ' formula',(Ln/vn)*(np.exp(vn*xn/(cn*Ln))-1))
r2=solve_ivp(lambda tt,y:[vn*y[0]/Ln-cn],[0,1e4],[0.05],events=ev,rtol=1e-10,atol=1e-12)
print('(a) x0=5cm numeric',r2.t_events[0],' formula',(Ln/vn)*np.log(2))
r3=solve_ivp(lambda tt,y:[vn*y[0]/Ln-cn],[0,500],[0.5],rtol=1e-10); print('(a) x0=50cm x(500s)=',r3.y[0,-1])
# (b) maximum distance from wall
r4=solve_ivp(lambda tt,y:[vn*y[0]/(Ln+vn*tt)-cn],[0,14741],[xn],rtol=1e-11,atol=1e-12,dense_output=True)
ts=np.linspace(0,14741,2000001); xs=r4.sol(ts)[0]
print('(b) max x numeric',xs.max(),'at t',ts[xs.argmax()],' formula',cn*Ln/vn*np.exp(vn*xn/(cn*Ln)-1),' t formula',(Ln/vn)*(np.exp(vn*xn/(cn*Ln)-1)-1))
# (4)(1): x0=50cm -> reaches winch x=L ; time
ev2=lambda tt,y:y[0]-Ln; ev2.terminal=True
r5=solve_ivp(lambda tt,y:[vn*y[0]/Ln-cn],[0,1e4],[xn],events=ev2,rtol=1e-11); print('(a) x0=50cm reaches winch at',r5.t_events[0],' formula',Ln/vn*np.log((Ln-0.1)/(xn-0.1)))
