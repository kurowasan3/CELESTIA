# Stacked-ball tower: sequential collisions, check mass ratios and heights
import sympy as sp
mu,e,v0=sp.symbols('mu e v0',positive=True)
m1=1; m2=mu
def elastic(ma,ua,mb,ub):
    # 1D elastic collision, returns (ua',ub')
    ua2=((ma-mb)*ua+2*mb*ub)/(ma+mb); ub2=((mb-ma)*ub+2*ma*ua)/(ma+mb)
    return sp.simplify(ua2),sp.simplify(ub2)
# (1) e=1, two balls, up positive
a,b=elastic(1,v0,mu,-v0)
print('(1) v1\'=',sp.factor(a),' v2\'=',sp.factor(b))
print('    H2/h=',sp.factor((b/v0)**2), ' v1\'=0 at mu=',sp.solve(sp.Eq(a,0),mu))
# (2) floor e
a,b=elastic(1,e*v0,mu,-v0)
sol=sp.solve(sp.Eq(a,0),mu); print('(2) mu=',sol)
b2=sp.simplify(b.subs(mu,sol[0])); print('    v2\'=',sp.factor(b2),' H2/h=',sp.factor((b2/v0)**2))
# (3) n-ball tower e=1: each lower ball stops
k=sp.symbols('k',positive=True,integer=True)
u=[None,sp.Integer(1)]   # in units of v0
ms=[None,sp.Integer(1)]
for i in range(1,8):
    r=sp.symbols('r',positive=True)
    a,bb=elastic(ms[i],u[i],ms[i]*r,-1)
    rs=[s for s in sp.solve(sp.Eq(a,0),r) if s>0][0]
    ms.append(ms[i]*rs); u.append(sp.simplify(bb.subs(r,rs)))
print('(3) u_k=',u[1:]); print('    m_k=',ms[1:])
print('    formula 2/(k(k+1)) ok:',all(ms[i]==sp.Rational(2,i*(i+1)) for i in range(1,9)))
# energy check: sum m_k * v0^2 == m_n * (n v0)^2
for n in range(1,9):
    assert sum(ms[1:n+1])==ms[n]*u[n]**2
print('    energy check ok; n=4: H/h=',u[4]**2,' m4/M=',ms[4]/sum(ms[1:5]))
# numeric event-driven simulation for n=4 with finite tiny gaps (independent check)
import numpy as np
g=9.8;h=1.0;eps=1e-4
mass=np.array([float(x) for x in ms[1:5]])
y=np.array([h+2*eps*i for i in range(4)]);v=np.zeros(4)
t=0;dt=2e-7
# simple integration with collision handling (point balls, gaps eps)
maxh=0
for step in range(int(0.5/dt)):
    v-=g*dt; y+=v*dt
    if y[0]<0 and v[0]<0: v[0]=-v[0]; y[0]=-y[0]
    for i in range(3):
        if y[i+1]-y[i]<0 and v[i+1]-v[i]<0:
            ma,mb=mass[i],mass[i+1];ua,ub=v[i],v[i+1]
            v[i]=((ma-mb)*ua+2*mb*ub)/(ma+mb); v[i+1]=((mb-ma)*ub+2*ma*ua)/(ma+mb)
print('    sim: velocities after collisions/v0 =',np.round(v/np.sqrt(2*g*h),3),' top height/h ~',round((y[3]+v[3]**2/(2*g))/h,3),' (expected 16)')
