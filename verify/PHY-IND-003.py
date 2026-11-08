# N38-1 verification: rod launched up inclined rails in vertical B with friction
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,g,B,l,R,th,mu,v,v0=sp.symbols('m g B l R theta mu v v0',positive=True)
F=B**2*l**2*v*sp.cos(th)/R          # horizontal magnetic force magnitude
# ascending: normal force, along-slope accel (up positive)
N_up=m*g*sp.cos(th)-F*sp.sin(th)
a_up=sp.simplify((-m*g*sp.sin(th)-F*sp.cos(th)-mu*N_up)/m)
print('N_up =',sp.factor(N_up)); print('vc =',sp.solve(N_up,v))
print('a_up =',sp.collect(sp.expand(a_up),v))
# special mu = cot(theta)
a_sp=sp.simplify(a_up.subs(mu,sp.cot(th)))
print('a_up(mu=cot) =',a_sp)
t1=v0/(g/sp.sin(th)); d=v0**2*sp.sin(th)/(2*g)
t=sp.symbols('t',positive=True); vt=v0-g/sp.sin(th)*t
Q=sp.integrate((B*l*vt*sp.cos(th))**2/R,(t,0,t1))
print('Q =',sp.simplify(Q))
# energy check: 1/2 m v0^2 = m g d sin + friction work + Q
Wf=sp.integrate(sp.cot(th)*(m*g*sp.cos(th)-B**2*l**2*vt*sp.cos(th)*sp.sin(th)/R)*vt,(t,0,t1))
print('energy residual =',sp.simplify(m*v0**2/2-m*g*d*sp.sin(th)-Wf-Q))
# descending terminal speed
u=sp.symbols('u',positive=True); Fu=B**2*l**2*u*sp.cos(th)/R
ug=sp.solve(sp.Eq(m*g*sp.sin(th),Fu*sp.cos(th)+mu*(m*g*sp.cos(th)+Fu*sp.sin(th))),u)[0]
us=sp.simplify(ug.subs(mu,sp.cot(th)))
print('u general =',sp.simplify(ug)); print('u(mu=cot) =',us)
print('check form:',sp.simplify(us-(-m*g*R*sp.cos(2*th)/(2*B**2*l**2*sp.sin(th)*sp.cos(th)**2))))
# numerical full process
vals=dict(m=0.2,g=9.8,B=0.8,l=0.5,R=0.5,th=np.radians(60))
mv,gv,Bv,lv,Rv,tv=vals['m'],vals['g'],vals['B'],vals['l'],vals['R'],vals['th']
muv=1/np.tan(tv); vcv=mv*gv*Rv/(Bv**2*lv**2*np.sin(tv)); v0v=0.8*vcv
def rhs(t,y):
    s,vv=y; Fh=Bv**2*lv**2*abs(vv)*np.cos(tv)/Rv
    if vv>0:
        N=mv*gv*np.cos(tv)-Fh*np.sin(tv); assert N>-1e-9
        a=(-mv*gv*np.sin(tv)-Fh*np.cos(tv)-muv*N)/mv
    else:
        N=mv*gv*np.cos(tv)+Fh*np.sin(tv)
        a=(-mv*gv*np.sin(tv)+Fh*np.cos(tv)+muv*N)/mv
    return [vv,a]
ev=lambda t,y:y[1]; ev.terminal=True; ev.direction=-1
s1=solve_ivp(rhs,[0,10],[0,v0v],events=ev,rtol=1e-10,atol=1e-12,max_step=1e-3)
print('num t1,d:',s1.t[-1],s1.y[0,-1],' formula:',v0v*np.sin(tv)/gv,v0v**2*np.sin(tv)/(2*gv))
s2=solve_ivp(rhs,[0,60],[0,-1e-9],rtol=1e-10,atol=1e-12)
print('num terminal:',s2.y[1,-1],' formula:',-(-mv*gv*Rv*np.cos(2*tv)/(2*Bv**2*lv**2*np.sin(tv)*np.cos(tv)**2)))
