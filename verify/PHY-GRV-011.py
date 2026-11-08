# J049-1: launch from Earth's surface with v0^2=gR at elevation alpha
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
x,al,g,R=sp.symbols('x alpha g R',positive=True)
# apogee r=xR: angular momentum R v0 cos(al) = xR va ; energy
eq=sp.Eq(sp.Rational(1,2)*g*R*sp.cos(al)**2/x**2 - g*R/x, sp.Rational(1,2)*g*R - g*R)
print('apogee x:',sp.solve(eq,x))
# semimajor axis from energy: v0^2/2 - gR = -gR^2/(2a)
a=sp.symbols('a',positive=True); print('a =',sp.solve(sp.Eq(g*R/2-g*R,-g*R**2/(2*a)),a))
# time: T*(1/2 + e/pi), e=sin(al), T=2pi sqrt(R/g)
t_formula=lambda alpha,g_,R_: np.sqrt(R_/g_)*(np.pi+2*np.sin(alpha))
g_=9.8;R_=6.4e6;GM=g_*R_**2
for alpha in [np.pi/4, 0.3, 1.2]:
    v0=np.sqrt(g_*R_)
    Y0=[0,R_, v0*np.cos(alpha), v0*np.sin(alpha)]  # launch at top point (0,R), horizontal = +x
    def f(t,Y):
        r3=np.hypot(Y[0],Y[1])**3; return [Y[2],Y[3],-GM*Y[0]/r3,-GM*Y[1]/r3]
    def land(t,Y): return np.hypot(Y[0],Y[1])-R_
    land.terminal=True; land.direction=-1
    s=solve_ivp(f,[1,1e5],[Y0[0]+Y0[2],Y0[1]+Y0[3],Y0[2],Y0[3]] if False else Y0,events=land,rtol=1e-11,atol=1e-6,first_step=1.0)
    # avoid immediate event at t=0: start integration from tiny time
    s=solve_ivp(f,[0,1e5],Y0,events=land,rtol=1e-11,atol=1e-6)
    te=s.t_events[0]; te=te[te>1]; ye=s.y_events[0][s.t_events[0]>1]
    tl=te[0]; xl,yl=ye[0][:2]
    ang=np.arccos((yl)/R_)
    rmax=np.max(np.hypot(s.y[0],s.y[1]))
    print(f'alpha={alpha:.3f}: t num={tl:.2f} formula={t_formula(alpha,g_,R_):.2f}; angle num={ang:.6f} formula={np.pi-2*alpha:.6f}; rmax/R={rmax/R_:.6f} vs {1+np.sin(alpha):.6f}')
print('sqrt(R/g)=',np.sqrt(R_/g_),' t(pi/4)=',t_formula(np.pi/4,g_,R_))
