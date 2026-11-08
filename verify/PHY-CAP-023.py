# S0504-1 verification: spring-supported plate + parallel reservoir capacitor, total charge Q fixed
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
x,d0,k,C0,beta,Q=sp.symbols('x d0 k C0 beta Q',positive=True)
eS=C0*d0                    # eps0*S
Cp=beta*C0
QA=Q*eS/(eS+Cp*x)
F=sp.simplify(QA**2/(2*eS)) # attraction (toward closing)
print('F =',F)
net=k*(d0-x)-F              # force in +x (opening) direction
# critical point: net=0 and d(net)/dx=0
Q2=sp.solve(sp.Eq(net,0),Q)[0]**2
Q2=sp.simplify(Q2)
dQ2=sp.diff(Q2,x)
xs=sp.solve(sp.Eq(dQ2,0),x); print('stationary pts of Q^2(x):',xs)
xstar=[s for s in xs if sp.simplify(s-d0)!=0][0]
print('x* =',sp.simplify(xstar))
assert sp.simplify(xstar-d0*(2*beta-1)/(3*beta))==0
Qc=sp.sqrt(sp.simplify(Q2.subs(x,xstar)))
print('Qc =',sp.simplify(Qc))
Qc_expect=sp.Rational(2,3)*(1+beta)**sp.Rational(3,2)*d0*sp.sqrt(2*k*C0/(3*beta))
assert sp.simplify(Qc**2-Qc_expect**2)==0
# book case beta=1 -> (8/3) d sqrt(kC0/3)
assert sp.simplify(Qc_expect.subs(beta,1)-sp.Rational(8,3)*d0*sp.sqrt(k*C0/3))==0
# stability: d(net)/dx<0 at equilibrium  <=> x > x*
stab=sp.simplify(sp.diff(net,x).subs(Q,sp.sqrt(Q2)))
print('d(net)/dx at eq =',sp.factor(stab))
# contact charge for beta<=1/2
Qcont=sp.sqrt(sp.limit(Q2,x,0)); print('Q at contact =',sp.simplify(Qcont))
assert sp.simplify(Qcont**2-2*k*C0*d0**2)==0
# numeric: quasi-static ramp of Q with damping, beta=1 and beta=0.3
def run(bv,Qmax,T=4000.0):
    kv,d0v,C0v,m,g=1.0,1.0,1.0,1.0,0.6
    def rhs(t,y):
        xx,v=y; Qv=Qmax*min(t/(0.8*T),1.0)
        Fv=Qv**2*(C0v*d0v)/(2*(C0v*d0v+bv*C0v*xx)**2)
        return [v,(kv*(d0v-xx)-Fv-g*v)/m]
    ev=lambda t,y: y[0]-1e-3; ev.terminal=True
    s=solve_ivp(rhs,[0,T],[1.0,0.0],events=ev,max_step=0.5,rtol=1e-8)
    return s.y[0,-1], (s.t_events[0][0] if len(s.t_events[0]) else None)
for bv in [1.0,2.0]:
    qc=float(Qc_expect.subs({beta:bv,d0:1,k:1,C0:1})); xsv=(2*bv-1)/(3*bv)
    xb,tb=run(bv,0.99*qc); xa,ta=run(bv,1.01*qc)
    print(f'beta={bv}: Qc={qc:.4f} x*={xsv:.4f} | 0.99Qc -> final x={xb:.4f} contact={tb} | 1.01Qc -> contact at t={ta}')
    xx_=sp.symbols('xx_'); xe=sp.nsolve((1-xx_)*2*(1+bv*xx_)**2/1-(0.99*qc)**2,xx_,xsv+0.1)
    print('   expected stable equilibrium at 0.99Qc:',xe)
    assert tb is None and ta is not None and abs(xb-float(xe))<0.01 and xe>xsv
bv=0.3; qk=float(sp.sqrt(2))
xb,tb=run(bv,0.95*qk); print('beta=0.3, Q=0.95*sqrt(2kC0)d0 -> x=',xb, 'expected eq x from cubic')
# expected equilibrium
xx=sp.symbols('xx'); sol=sp.nsolve((1-xx)*2*(1+0.3*xx)**2-(0.95*qk)**2,xx,0.05); print('cubic eq x=',sol)
assert abs(xb-float(sol))<0.01
print('ALL OK')
