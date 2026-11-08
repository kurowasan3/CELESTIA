# J026-1 verification: block on quarter circle on cart moving at constant u
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,R,u,m,phi,H=sp.symbols('g R u m phi H',positive=True)
c=sp.symbols('c',positive=True)
# cart frame (inertial): angle phi from vertical below O; start phi=pi/2 at P
w2=2*g*R*sp.cos(phi)              # relative speed^2
N=m*g*sp.cos(phi)+m*w2/R           # normal force
print('N(phi)=',sp.simplify(N),' N_B=',N.subs(phi,0))
wB=sp.sqrt(2*g*R)
# ground frame energy theorem: block velocity at B = u - wB (rightward positive)
WN=sp.Rational(1,2)*m*(u-wB)**2-sp.Rational(1,2)*m*u**2 - m*g*R
print('W_N =',sp.simplify(WN))
# impulse method: horizontal impulse = m*(-wB); W=u*J
print('u*J =',sp.simplify(u*m*(-wB)), sp.simplify(WN-u*m*(-wB)))
# (4) u=sqrt(2gR): ground speed^2 along path, relative velocity direction (-cos, -sin)
v2=(u-sp.sqrt(2*g*R*c)*c)**2+(2*g*R*c)*(1-c**2)
v2s=sp.simplify(v2.subs(u,sp.sqrt(2*g*R)))
print('v^2(c)=',sp.expand(v2s))
crit=sp.solve(sp.diff(v2s,c),c); print('crit',crit, [sp.simplify(v2s.subs(c,cc)) for cc in crit])
print('landing offset', sp.simplify(sp.sqrt(2*g*R)*sp.sqrt(2*H/g)))
# numerical: integrate in cart frame, compute work of N in ground frame
gv,Rv,uv,mv=9.8,0.5,2.0,1.0
def f(t,y):
    ph,om=y; return [om, -gv/Rv*np.sin(ph)]
ev=lambda t,y:y[0]; ev.terminal=True; ev.direction=-1
sol=solve_ivp(f,[0,5],[np.pi/2,0],events=ev,max_step=1e-4,rtol=1e-10,atol=1e-12)
ph,om=sol.y; t=sol.t
# position rel O: (R sin ph, -R cos ph); velocity rel = R om (cos ph, sin ph)
vrx=Rv*om*np.cos(ph); vry=Rv*om*np.sin(ph)
Nmag=mv*gv*np.cos(ph)+mv*Rv*om**2
# normal direction on block: toward center O: (-sin ph, cos ph)
Nx=-Nmag*np.sin(ph); Ny=Nmag*np.cos(ph)
P=Nx*(uv+vrx)+Ny*vry
W=np.trapezoid(P,t); print('numeric W_N',W,' formula',-mv*uv*np.sqrt(2*gv*Rv))
# u = sqrt(2gR): max ground speed
uv=np.sqrt(2*gv*Rv); vg2=(uv+vrx)**2+vry**2
i=np.argmax(vg2); print('max ground v^2/(2gR)',vg2[i]/(2*gv*Rv),'cos at max',np.cos(ph[i]),' expected',28/27,1/9)
print('final ground speed', np.sqrt(vg2[-1]))
