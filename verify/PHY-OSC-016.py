import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
rho,S,g,L,th,s,P,l,r=sp.symbols('rho S g L theta s P l r',positive=True)
# potential energy method: move slug of length s from right top to left top
# heights: left surface at h=0 rises; right surface along incline. Energy of slug: centre from -s sin(th)/2 to +s/2
U=rho*S*s*g*(s/2+s*sp.sin(th)/2)
F=-sp.diff(U,s); print('F',sp.simplify(F))
w0=sp.sqrt(sp.diff(U,s,2)/(rho*S*L)); tau0=2*sp.pi/w0; print('tau0',sp.simplify(tau0))
# closed arm isothermal: P l = (P+dP)(l-s)
dP=sp.series(P*l/(l-s)-P,s,0,2).removeO(); print('dP',dP)
k1=rho*S*g*(1+sp.sin(th))+S*sp.diff(dP,s); tau1=2*sp.pi*sp.sqrt(rho*S*L/k1); print('tau1',sp.simplify(tau1))
Psol=sp.solve(sp.Eq((tau0/tau1)**2,r**2),P); print('P',Psol)
val=Psol[0].subs({rho:1000,g:9.8,l:0.20,th:sp.pi/6,r:5.9}); print('P num',sp.N(val))
# numeric ODE nonlinear isothermal, small amplitude
vals=dict(rho=1000,g=9.8,L=0.8,th=np.pi/6,P=1.0e5,l=0.2,S=1e-4)
def f(t,y):
    s_,v=y; d=vals
    Fg=-d['rho']*d['S']*d['g']*(1+np.sin(d['th']))*s_
    Fp=-(d['P']*d['l']/(d['l']-s_)-d['P'])*d['S']
    return [v,(Fg+Fp)/(d['rho']*d['S']*d['L'])]
ev=lambda t,y:y[1]; ev.direction=-1
sol=solve_ivp(f,[1e-7,5],[1e-5,-1e-12],events=ev,rtol=1e-11,atol=1e-16)
Tnum=sol.t_events[0][0]
Tan=float(tau1.subs({rho:1000,g:9.8,L:0.8,th:sp.pi/6,P:1e5,l:0.2}))
print('tau1 numeric',Tnum,'analytic',Tan); assert abs(Tnum-Tan)/Tan<1e-4; print('OK')
