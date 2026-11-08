# S0528-1: rod on rails in loop encircling solenoid; flux Phi_s changes 0->Phi1
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
# equations: R I = -dPhis/dt - B l u ; m du/dt = B l I  (I counterclockwise from above, x to the right)
m,R,B,l,P1,T,t=sp.symbols('m R B l Phi1 T t',positive=True)
tm=m*R/(B*l)**2; a=P1/T
u=sp.Function('u')
us=sp.dsolve(sp.Eq(m*u(t).diff(t),B*l*(-a-B*l*u(t))/R),u(t),ics={u(0):0}).rhs
uexp=-(a/(B*l))*(1-sp.exp(-t/tm))
print('(5) u(t)',sp.simplify(us-uexp))
uT=uexp.subs(t,T)
x_ramp=sp.integrate(uexp,(t,0,T)); x_after=uT*tm
print('(5) total displacement',sp.simplify(x_ramp+x_after))  # expect -Phi1/(B l)
I_ramp=sp.simplify(m*uexp.diff(t)/(B*l)); print('I during ramp',I_ramp)
# fast limit: impulse
v1=B*l*P1/(m*R); print('(fast) disp',sp.simplify(v1*tm), ' KE',sp.simplify(m*v1**2/2))
# generic profiles numeric
def run(prof, dprof, mv=0.2,Rv=0.5,Bv=0.8,lv=0.3,P1v=0.06):
    def f(tt,y):
        x,uu,q=y
        I=(-dprof(tt)-Bv*lv*uu)/Rv
        return [uu,Bv*lv*I/mv,I]
    s=solve_ivp(f,[0,60],[0,0,0],rtol=1e-10,atol=1e-13,max_step=0.005)
    return s.y[0,-1],s.y[1,-1],s.y[2,-1],-P1v/(Bv*lv)
P1v=0.06
for name,dp in [('exp tau=0.3',lambda tt:P1v/0.3*np.exp(-tt/0.3)),
                ('sin^2 pulse 2s',lambda tt:(P1v/2*(1-np.cos(np.pi*tt/2))*0+ (P1v*np.pi/4*np.sin(np.pi*tt/2) if tt<2 else 0))),
                ('linear 1s',lambda tt:P1v/1.0 if tt<1 else 0.0)]:
    print(name,'x_final,u_final,Q =',run(None,dp),' expect x=-Phi1/(Bl)')
# sign of current change exactly at T for linear ramp (numeric)
mv,Rv,Bv,lv,Tv=0.2,0.5,0.8,0.3,1.0
def f(tt,y):
    dP=P1v/Tv if tt<Tv else 0.0
    I=(-dP-Bv*lv*y[1])/Rv; return [y[1],Bv*lv*I/mv]
s=solve_ivp(f,[0,5],[0,0],rtol=1e-10,atol=1e-13,max_step=0.001,dense_output=True)
for tt in [0.5,0.999,1.001,2]:
    y=s.sol(tt); dP=P1v/Tv if tt<Tv else 0.0
    print('t',tt,'I',(-dP-Bv*lv*y[1])/Rv)
