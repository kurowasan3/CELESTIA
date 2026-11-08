# K20012-1 verification: fountain with absorbing/reflecting plate, arrival-time edge
import sympy as sp, numpy as np
V,g,D=sp.symbols('V g D',positive=True)
VD=sp.sqrt(2*g*D)
u=sp.sqrt(V**2-2*g*D)
t_ref=2*(V-u)/g
# check by explicit kinematics: rise to D time t1, fall from D with speed u: D = u t2 + g t2^2/2
t1=(V-u)/g
t2=sp.solve(sp.Eq(D,u*sp.Symbol('t2')+g*sp.Symbol('t2')**2/2),sp.Symbol('t2'))
t2=[x for x in t2 if sp.simplify(x.subs({V:3,g:1,D:1}))>0][0]
print('t_ref check:',sp.simplify(t1+t2-t_ref))
print('d t_ref/dV =',sp.simplify(sp.diff(t_ref,V)))
# V for t = t_c/2
tc=2*VD/g
print('V at t_c/2:',sp.solve(sp.Eq(t_ref,tc/2),V))
# numeric: g from given D,tc
Dn=0.800; tcn=0.8080
gn=8*Dn/tcn**2; print('g =',gn)
# Monte Carlo: uniform V in [0,5/4 VD], fraction arriving in [tc/2, tc]
gg=9.8; DD=0.8; vd=np.sqrt(2*gg*DD); tcc=2*vd/gg
Vs=np.random.default_rng(1).uniform(0,1.25*vd,2_000_000)
ts=np.where(Vs<vd,2*Vs/gg,2*(Vs-np.sqrt(np.maximum(Vs**2-2*gg*DD,0)))/gg)
print('fraction in [tc/2,tc]:',np.mean((ts>=tcc/2)&(ts<=tcc)),' expected 0.6')
print('max t / tc:',ts.max()/tcc,'  min t of reflected / tc:',ts[Vs>vd].min()/tcc)
# fast simulation check of reflected time with event-based integration for one V
from scipy.integrate import solve_ivp
Vt=1.2*vd
def rhs(t,y): return [y[1],-gg]
ev=lambda t,y: y[0]-DD; ev.terminal=True; ev.direction=1
s1=solve_ivp(rhs,[0,5],[0,Vt],events=ev,rtol=1e-12,atol=1e-12,max_step=1e-3)
th=s1.t_events[0][0]; vh=s1.y_events[0][0][1]
ev2=lambda t,y: y[0]; ev2.terminal=True; ev2.direction=-1
s2=solve_ivp(rhs,[0,5],[DD,-vh],events=ev2,rtol=1e-12,atol=1e-12,max_step=1e-3)
print('ivp reflected time',th+s2.t_events[0][0],' formula',float(t_ref.subs({V:Vt,g:gg,D:DD})))
