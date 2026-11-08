import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
p,V,dV,Cv,R=sp.symbols('p V dV C_V R',positive=True)
# first law adiabatic: Cv n dT = -p dV, state: p dV + V dp = nR dT
dp,ndT=sp.symbols('dp ndT')
sol=sp.solve([sp.Eq(Cv*ndT,-p*dV),sp.Eq(p*dV+V*dp,R*ndT)],[dp,ndT])
print('dp',sp.simplify(sol[dp]))
P,S,M,V0,V1,x=sp.symbols('P S M V0 V1 x',positive=True)
g1=sp.Rational(5,3); g2=sp.Rational(7,5)
F=(-g1*P*S*x/V1 - g2*P*S*x/(V0-V1))*S
K=sp.simplify(-F/x); tau=2*sp.pi*sp.sqrt(M/K); print('tau',sp.simplify(tau))
crit=sp.solve(sp.diff(g1/V1+g2/(V0-V1),V1),V1); print('crit',crit)
V1s=5*V0/(5+sp.sqrt(21)); print('check', [sp.simplify(c-V1s) for c in crit])
taumax=sp.simplify(tau.subs(V1,V1s)); print('taumax',taumax, sp.simplify(taumax/(2*sp.pi*sp.sqrt(15)/(5+sp.sqrt(21))/S*sp.sqrt(M*V0/P))))
# numeric: nonlinear adiabatic, units P=S=M=V0=1
def period(v1,amp=1e-4):
    v2=1-v1
    def f(t,y):
        X,vv=y; pl=(v1/(v1+X))**(5/3); pr=(v2/(v2-X))**(7/5); return [vv,(pl-pr)]
    ev=lambda t,y:y[1]; ev.direction=1
    s=solve_ivp(f,[0,50],[amp,0],events=ev,rtol=1e-11,atol=1e-14)
    return s.t_events[0][0]  # starts at max displacement, velocity crosses 0 upward at the min? -> half? 
def period2(v1,amp=1e-4):
    v2=1-v1
    def f(t,y):
        X,vv=y; pl=(v1/(v1+X))**(5/3); pr=(v2/(v2-X))**(7/5); return [vv,(pl-pr)]
    ev=lambda t,y:y[1]; ev.direction=-1
    s=solve_ivp(f,[1e-6,50],[amp,-1e-14],events=ev,rtol=1e-11,atol=1e-16)
    return s.t_events[0][0]
for v1 in [0.3,0.5,float(5/(5+np.sqrt(21))),0.6]:
    Tn=period2(v1); Ta=float(tau.subs({P:1,S:1,M:1,V0:1,V1:v1}))
    print(v1,Tn,Ta)
vs=np.linspace(0.3,0.7,401); Ts=[float(tau.subs({P:1,S:1,M:1,V0:1,V1:v})) for v in vs]
print('argmax',vs[int(np.argmax(Ts))], float(5/(5+np.sqrt(21))))
assert abs(period2(0.5)-float(tau.subs({P:1,S:1,M:1,V0:1,V1:0.5})))<1e-5
assert abs(vs[int(np.argmax(Ts))]-0.5218)<1e-3; print('OK')
