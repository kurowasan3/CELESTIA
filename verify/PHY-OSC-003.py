# M28-1: string length l from A, released from horizontal; peg B at depth x below A; after catching, radius R=l-x.
# Find x such that after slack the ball hits the peg.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
l,g,R,s,m=sp.symbols('l g R s m',positive=True)
# slack point at angle alpha above horizontal through peg: v^2 = g R s, energy from release (A level)
eqE=sp.Eq(g*R*s, 2*g*(l-R-R*s))
s_of_R=sp.solve(eqE,s)[0]; print('sin alpha =',s_of_R)
# projectile through peg -> s=1/sqrt3 (derived independently)
c=sp.sqrt(1-s**2); t=sp.symbols('t',positive=True); v=sp.sqrt(g*R*s)
X=R*c - v*s*t; Y=R*s + v*c*t - g*t**2/2
tt=sp.solve(sp.Eq(X,0),t)[0]
ssol=sp.solve(sp.Eq(sp.simplify(Y.subs(t,tt)),0),s); print('s for hitting peg:',ssol)
Rsol=sp.solve(sp.Eq(s_of_R,1/sp.sqrt(3)),R)[0]; print('R=',sp.nsimplify(sp.radsimp(Rsol)),' x=l-R=',sp.radsimp(l-Rsol))
xsol=sp.radsimp(l-Rsol)
assert sp.simplify(xsol-(2*sp.sqrt(3)-3)*l)==0
# slack occurs above peg level iff 0<s<1 -> range of x
print('range: s<1 ->', sp.solve(s_of_R<1,R), ' s>0 -> R<l')
# tension before/after catching
Tb=m*g+m*2*g*l/l; Ta=sp.radsimp(m*g+m*2*g*l/Rsol); print('T before',Tb,' after',sp.simplify(Ta))
# numeric simulation: l=1, g=9.8
lv,gv=1.0,9.8; Rv=float(Rsol.subs(l,1)); 
# phase about peg: angle phi from downward vertical, ccw; phi'' = -(g/R) sin phi; tension T/m = R phi'^2 + g cos phi
v0=np.sqrt(2*gv*lv)
def f(t,y): return [y[1], -(gv/Rv)*np.sin(y[0])]
def slack(t,y): return Rv*y[1]**2+gv*np.cos(y[0])
slack.terminal=True; slack.direction=-1
sol=solve_ivp(f,[0,10],[0,v0/Rv],events=slack,rtol=1e-12,atol=1e-12)
ph,om=sol.y_events[0][0]
print('slack angle above horizontal: sin=',np.sin(ph-np.pi/2),' 1/sqrt3=',1/np.sqrt(3))
# position rel peg: (R sin ph, -R cos ph), velocity R om (cos ph, sin ph)
p=np.array([Rv*np.sin(ph),-Rv*np.cos(ph)]); vv=Rv*om*np.array([np.cos(ph),np.sin(ph)])
def fp(t,y): return [y[2],y[3],0,-gv]
sol2=solve_ivp(fp,[0,2],[*p,*vv],rtol=1e-12,atol=1e-12,dense_output=True)
ts=np.linspace(0,2,200001); Y=sol2.sol(ts); dist=np.hypot(Y[0],Y[1])
i=np.argmin(dist); print('min distance to peg=',dist[i],' at t=',ts[i]); assert dist[i]<1e-4
assert np.all(dist[:i]<=Rv+1e-9)   # string stays slack (ball inside circle) until hitting peg
print('OK')
