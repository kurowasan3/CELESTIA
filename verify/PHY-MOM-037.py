# S0108-1: body on rough floor dropping cargo that ends at rest w.r.t. floor
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m0,v0,mu,g,t,dm,m,v = sp.symbols('m0 v0 mu g t dm m v', positive=True)
dv = sp.symbols('dv')
# (1) momentum of remainder: (m-dm)(v+dv) + dm*0 - m v = -mu m g dt  (to first order)
dt = sp.symbols('dt', positive=True)
sol = sp.solve(sp.Eq(sp.expand((m-dm)*(v+dv) - m*v).subs(dm*dv,0), -mu*m*g*dt), dv)[0]
print('dv =', sp.simplify(sol))
# Strategy A: constant speed -> dv=0 -> dm = mu m g dt / v0
mf = sp.Function('mm')
ode = sp.Eq(mf(t).diff(t), -mu*g*mf(t)/v0)
mA = sp.dsolve(ode, ics={mf(0):m0}).rhs; print('m_A(t)=',mA)
t1 = sp.solve(sp.Eq(mA, m0/2), t)[0]; print('t1=',t1)
LA = sp.simplify(v0*t1 + v0**2/(2*mu*g)); print('L_A =', LA)
WA = sp.simplify((sp.Rational(1,2)*(m0/2)*v0**2 - sp.Rational(1,2)*m0*v0**2) + sp.integrate(mu*g*mA*v0,(t,0,t1)))
print('W_A =', WA)
# Strategy B: eject m0/2 at once
vB = 2*v0; LB = vB**2/(2*mu*g); WB = sp.Rational(1,2)*(m0/2)*vB**2 - sp.Rational(1,2)*m0*v0**2
print('L_B =', LB, ' W_B =', WB)
# numeric: general strategy (eject at constant rate r until m0/2) vs B, A
gv,muv,v0v,m0v=9.8,0.2,5.0,1.0
def run(rate_fn):
    # state x, p, m ; dp/dt = -mu m g ; dm/dt = -rate (while m>m0/2 and moving)
    def f(t,y):
        x,p,mm=y
        vv=p/mm
        r = rate_fn(t,mm) if mm>m0v/2+1e-12 else 0.0
        return [vv, -muv*mm*gv, -r]
    stop=lambda t,y: y[1]; stop.terminal=True; stop.direction=-1
    s=solve_ivp(f,[0,100],[0,m0v*v0v,m0v],events=stop,rtol=1e-10,atol=1e-12,max_step=1e-3)
    return s.y_events[0][0][0]
LAn = run(lambda t,mm: muv*gv*mm/v0v)
print('numeric L_A:',LAn,' formula:',float(LA.subs({v0:v0v,mu:muv,g:gv})))
print('formula L_B:',float(LB.subs({v0:v0v,mu:muv,g:gv})))
for r in [0.5,2,10,100]:
    print(' const-rate',r,'-> L =',run(lambda t,mm,r=r: r))
# energy of device for constant-speed strategy numerically
