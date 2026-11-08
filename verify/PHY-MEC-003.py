import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
S,s,g,H,Q,h,t=sp.symbols('S s g H Q h t',positive=True)
T1=S*H/Q
# drain only
u=sp.Function('u')
T2=sp.simplify(2*sp.sqrt(H)/(s*sp.sqrt(2*g)/S)); print('T2',T2)
hinf=sp.solve(sp.Eq(Q,s*sp.sqrt(2*g*h)),h)[0]
print('hinf/H in T1,T2:',sp.simplify(hinf/H - (T2/(2*T1))**2))
# numeric with T1=3, T2=2 (minutes), H=1
H0=1.0; T1v=3.0; T2v=2.0
alpha=H0/T1v; beta=2*np.sqrt(H0)/T2v
f=lambda tt,y:[alpha-beta*np.sqrt(max(y[0],0))]
ev=lambda tt,y:y[0]-0.25; ev.terminal=True
sol=solve_ivp(f,[0,20],[1.0],events=ev,rtol=1e-11,atol=1e-13)
print('time full->H/4 both open',sol.t_events[0][0],'formula',(T2v**2/(2*T1v))*(1.5+2*np.log(2)))
x0,x1=3.0,1.5
F=lambda x:-x-np.log(abs(1-x))
print('formula2',(2*alpha/beta**2)*(F(x1)-F(x0)),'exact expr',(T2v**2/(2*T1v))*(1.5+2*np.log(2)))
sol2=solve_ivp(f,[0,200],[0.0],rtol=1e-11,atol=1e-13); print('h_inf numeric',sol2.y[0][-1],'pred',(T2v/(2*T1v))**2)
# drain only to H/4
g2=lambda tt,y:[-beta*np.sqrt(max(y[0],0))]
sol3=solve_ivp(g2,[0,20],[1.0],events=ev,rtol=1e-11); print('drain only to H/4',sol3.t_events[0][0])
