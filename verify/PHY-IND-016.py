# T20032-1 crossed-rail two-bar induction: symbolic + numerical (staged) check
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,R,B,l,v0,t=sp.symbols('m R B l v0 t',positive=True)
# stage 1 (crossed): I = Bl(u+v)/(2R); both forces -BlI (same direction)
u=sp.Function('u'); v=sp.Function('v')
I=B*l*(u(t)+v(t))/(2*R)
sol=sp.dsolve([sp.Eq(m*u(t).diff(t),-B*l*I),sp.Eq(m*v(t).diff(t),-B*l*I)],ics={u(0):v0,v(0):0})
us=sp.simplify(sol[0].rhs); vs=sp.simplify(sol[1].rhs)
print('u(t)=',us,' v(t)=',vs)
uinf=sp.limit(us,t,sp.oo); vinf=sp.limit(vs,t,sp.oo)
print('u_inf,v_inf =',uinf,vinf); assert sp.simplify(uinf-v0/2)==0 and sp.simplify(vinf+v0/2)==0
print('u-v conserved:',sp.simplify(us-vs))
I0=I.subs({u(t):v0,v(t):0}); print('I0=',I0)
Q1=sp.simplify(m*v0**2/2-m*uinf**2/2-m*vinf**2/2); print('Q1=',Q1)
# joule heat by integration
Qint=sp.integrate(2*R*(B*l*(us+vs)/(2*R))**2,(t,0,sp.oo)); print('Q1(int)=',sp.simplify(Qint))
# stage 3 (same section)
w=sp.symbols('w')  # relative velocity
tau=m*R/(B**2*l**2)
rel=sp.integrate(v0*sp.exp(-t/tau),(t,0,sp.oo)); print('approach distance=',sp.simplify(rel))
# numeric full simulation with positions
mv,Rv,Bv,lv,V0=1.0,1.0,1.0,1.0,1.0
gap0,gapw=0.0,0.3   # insulator from x=0 to 0.3 ; left section x<0, right x>0.3
x1_0=-30.0; L=1.5  # bar2 position chosen so that when bar1 leaves insulator distance ~ L
def section(x):
    return 'L' if x<gap0 else ('I' if x<gap0+gapw else 'R')
def rhs(tt,y):
    x1,u_,x2,v_=y
    s1,s2=section(x1),section(x2)
    if 'I' in (s1,s2): F1=F2=0.0
    elif s1==s2:
        Ic=Bv*lv*(u_-v_)/(2*Rv); F1=-Bv*lv*Ic; F2=+Bv*lv*Ic
    else:
        Ic=Bv*lv*(u_+v_)/(2*Rv); F1=-Bv*lv*Ic; F2=-Bv*lv*Ic
    return [u_,F1/mv,v_,F2/mv]
# choose x2_0 so that at moment bar1 exits insulator, bar2 is L ahead: approx (bar1 approach rate v0 constant until then)
# in stage1+2 relative velocity = v0 exactly, so separation decreases by v0*t
x2_0=40.0
ev=lambda tt,y: y[0]-(gap0+gapw); ev.terminal=True; ev.direction=1
s=solve_ivp(rhs,[0,500],[x1_0,V0,x2_0,0.0],events=ev,max_step=0.01,rtol=1e-10,atol=1e-12)
x1,u1,x2,v2=s.y[:,-1]; print('at exit of insulator: u,v=',u1,v2,' separation=',x2-x1)
sep_exit=x2-x1
s2=solve_ivp(rhs,[s.t[-1],s.t[-1]+60],s.y[:,-1],max_step=0.01,rtol=1e-10,atol=1e-12)
x1f,uf,x2f,vf=s2.y[:,-1]; print('final u,v=',uf,vf,' final sep=',x2f-x1f,' expected',sep_exit-V0*mv*Rv/(Bv**2*lv**2))
assert abs(uf)<1e-6 and abs(vf)<1e-6 and abs((x2f-x1f)-(sep_exit-1.0))<1e-3
# total momentum change of bars = -m v0
print('impulse on bars total =',mv*(uf+vf)-mv*V0)
print('OK')
