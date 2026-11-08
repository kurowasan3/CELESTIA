# M01-1 verification: projectile P vs block Q sliding on smooth incline
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,D,v,al,ph,t=sp.symbols('g D v alpha phi t',positive=True)
# ground frame
eS=sp.Matrix([sp.cos(ph),sp.sin(ph)]); eN=sp.Matrix([-sp.sin(ph),sp.cos(ph)])
P=v*t*sp.Matrix([sp.cos(ph+al),sp.sin(ph+al)])+sp.Matrix([0,-g*t**2/2])
Q=D*eS - sp.Rational(1,2)*g*sp.sin(ph)*t**2*eS
rel=sp.simplify(P-Q)
# (1) relative acceleration
arel=sp.simplify(sp.diff(rel,t,2)); print('a_rel',arel.T, 'dot eS',sp.simplify(arel.dot(eS)),'dot eN',sp.simplify(arel.dot(eN)))
# (2) normal comp zero at t1, tangential zero at t1
n=sp.simplify(rel.dot(eN)); s=sp.simplify(rel.dot(eS))
t1=sp.solve(sp.Eq(n/t,0),t)[0]; print('t1',sp.simplify(t1))
cond=sp.simplify(s.subs(t,t1)); print('tangential at t1 (=0 needed):',sp.simplify(sp.trigsimp(cond)))
print('check v^2 sin2a = gDcos: ', sp.simplify(cond.subs(D, v**2*sp.sin(2*al)/(g*sp.cos(ph)))))
# (3) fixed target min speed: range up the slope
R=sp.simplify((P.dot(eS)).subs(t, sp.solve(sp.Eq(P.dot(eN)/t,0),t)[0]))
print('range fixed:',sp.simplify(R))
Rmax=sp.simplify(R.subs(al,sp.pi/4-ph/2)); print('Rmax',sp.simplify(Rmax), ' expected v^2/(g(1+sin))', sp.simplify(Rmax - v**2/(g*(1+sp.sin(ph)))))
# numeric check of max over alpha
import math
for phv in [0.3,0.6]:
    f=sp.lambdify(al,R.subs({v:1,g:1,ph:phv}))
    A=np.linspace(0.001,math.pi/2-phv-0.001,20001); print(phv,max(f(A)),1/(1+math.sin(phv)))
# (4) collision point at v1, alpha=45
v1=sp.sqrt(g*D*sp.cos(ph)); tc=t1.subs({v:v1,al:sp.pi/4})
pos=sp.simplify((Q.dot(eS)).subs(t,tc)); print('collision pt along slope',sp.simplify(pos), sp.simplify(pos-D*(1-sp.tan(ph))))
# numeric ODE: integrate both, find min distance and check perpendicular height >=0
for phv,Dv in [(0.4,1.0),(1.0,1.0)]:
    gv=9.8; v1v=math.sqrt(gv*Dv*math.cos(phv)); a=math.pi/4
    def f(tt,y): return [y[2],y[3],0,-gv, y[6],y[7], -gv*math.sin(phv)*math.cos(phv), -gv*math.sin(phv)**2]
    y0=[0,0,v1v*math.cos(phv+a),v1v*math.sin(phv+a), Dv*math.cos(phv),Dv*math.sin(phv),0,0]
    T=2*v1v*math.sin(a)/(gv*math.cos(phv))
    sol=solve_ivp(f,[0,T],y0,rtol=1e-11,atol=1e-12,dense_output=True)
    yy=sol.y[:,-1]; print('phi',phv,'gap at T',math.hypot(yy[0]-yy[4],yy[1]-yy[5]),'pos along slope',yy[0]*math.cos(phv)+yy[1]*math.sin(phv),'pred',Dv*(1-math.tan(phv)))
    ts=np.linspace(0,T,2001); Y=sol.sol(ts); h=-Y[0]*math.sin(phv)+Y[1]*math.cos(phv); print(' min height above slope',h.min(), 'launch angle from horiz deg',math.degrees(phv+a))
