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

# ---- (5) limits ----
gs,Rs,r,th=sp.symbols('g R r theta',positive=True)
tsym=sp.sqrt(Rs/gs)*(sp.pi+2*sp.sin(al))
print('(5a) alpha->pi/2:',sp.limit(tsym,al,sp.pi/2))
v=sp.sqrt(gs*Rs*(2*Rs/r-1))
# energy check
assert sp.simplify(v**2/2-gs*Rs**2/r-(gs*Rs/2-gs*Rs))==0
I1=sp.integrate(4*Rs*sp.sin(th)**2,(th,sp.pi/4,sp.pi/2)); print('(5a) substituted integral =',sp.simplify(I1))
assert sp.simplify(I1-Rs*(sp.pi/2+1))==0
from scipy.integrate import quad
tup=quad(lambda rr: 1/np.sqrt(g_*R_*(2*R_/rr-1)),R_,2*R_)[0]
print('(5a) numeric round trip=%.3f formula=%.3f'%(2*tup,(np.pi+2)*np.sqrt(R_/g_)))
assert abs(2*tup-(np.pi+2)*np.sqrt(R_/g_))<1e-3*2*tup
# radial launch simulation
def f(t,Y):
    r3=np.hypot(Y[0],Y[1])**3; return [Y[2],Y[3],-GM*Y[0]/r3,-GM*Y[1]/r3]
def land(t,Y): return np.hypot(Y[0],Y[1])-R_
land.terminal=True; land.direction=-1
s=solve_ivp(f,[0,1e5],[0,R_,0,np.sqrt(g_*R_)],events=land,rtol=1e-11,atol=1e-6)
print('(5a) radial sim landing time=%.3f'%s.t_events[0][-1])
print('(5b) alpha->0:',sp.limit(tsym,al,0),' angle->',sp.limit(sp.pi-2*al,al,0))
for alpha in [1e-2,1e-3]:
    v0=np.sqrt(g_*R_)
    s=solve_ivp(f,[0,1e5],[0,R_,v0*np.cos(alpha),v0*np.sin(alpha)],events=land,rtol=1e-12,atol=1e-6)
    te=s.t_events[0]; te=te[te>1]
    print('(5b) alpha=%g landing t/(pi sqrt(R/g))=%.5f'%(alpha,te[0]/(np.pi*np.sqrt(R_/g_))))
