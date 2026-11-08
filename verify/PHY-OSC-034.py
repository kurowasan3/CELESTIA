# J051-1: ball in smooth hemispherical bowl, launched horizontally (tangential) off-equilibrium
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
c,c0,w,g,R,m=sp.symbols('c c0 w g R m',positive=True)
s0sq=1-c0**2
# turning points: w*s0^2 = (1-c^2)(w + 2(c-c0)),  w=v0^2/(gR), c=cos(phi) depth=R c
print('circular w:', sp.solve(sp.diff((1-c**2)*(w+2*(c-c0))-w*s0sq,c).subs(c,c0),w))
print('reach rim (c=0):', sp.solve(sp.Eq(w*s0sq,(w+2*(0-c0))),w))
expr=sp.expand(w*s0sq-(1-c**2)*(w+2*(c-c0)))
for wv in [3]:
    e=expr.subs({c0:sp.Rational(1,2),w:wv}); print('w=',wv,sp.factor(e), sp.solve(e,c))
cs=(sp.sqrt(13)-3)/4
v2=3*g*R*sp.Rational(3,4)/(1-cs**2)
N=sp.nsimplify(sp.simplify(m*v2/R+m*g*cs)); print('v^2 top=',sp.radsimp(sp.simplify(v2)),' N=',sp.radsimp(N), float(N.subs({m:1,g:1})))
print('v2 numeric/gR', float(v2.subs({g:1,R:1})))
# numeric spherical pendulum: coords (phi from downward vertical, psi azimuth)
g_=9.8;R_=1.0
phi0=np.pi/3; v0=np.sqrt(3*g_*R_)
L=R_*np.sin(phi0)*v0  # = r_perp * v_psi  (per unit mass, /R for ang)
def f(t,y):
    ph,dph=y
    s=np.sin(ph); psid=L/(R_**2*s**2)
    return [dph, s*np.cos(ph)*psid**2 - g_/R_*s]
s=solve_ivp(f,[0,5],[phi0,0],rtol=1e-11,atol=1e-12,max_step=1e-3)
cmin=np.cos(s.y[0]).min(); print('numeric min cos phi (highest)=',cmin,' formula',(np.sqrt(13)-3)/4)
# N along path >0 check
ph=s.y[0]; v2n=(R_*s.y[1])**2+(L/(R_*np.sin(ph)))**2
Nn=v2n/R_+g_*np.cos(ph); print('min N/m =',Nn.min(),' at highest N/(mg) formula',float(N.subs({m:1,g:1})))
# rim case: c0=1/2, w=4
L2=R_*np.sin(phi0)*np.sqrt(4*g_*R_)
def f2(t,y):
    ph,dph=y; s_=np.sin(ph); psid=L2/(R_**2*s_**2); return [dph, s_*np.cos(ph)*psid**2-g_/R_*s_]
s2=solve_ivp(f2,[0,3],[phi0,0],rtol=1e-11,atol=1e-12,max_step=1e-3); print('rim case max phi(deg)=',np.degrees(s2.y[0].max()))
