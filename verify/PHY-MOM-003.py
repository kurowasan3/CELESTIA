# M14-1: projectile (2m) splits at apex into m + m; P1 lands at D (below apex) in half the no-split fall time.
import sympy as sp, numpy as np
v0,th,g,m,u,t = sp.symbols('v0 theta g m u t', positive=True)
V = v0*sp.cos(th); T = v0*sp.sin(th)/g; h = g*T**2/2
# P1 lands at D after T/2 => horizontal 0, vertical down speed u
u_s = sp.solve(sp.Eq(h - u*T/2 - g*(T/2)**2/2, 0), u)[0]; print('(1) u=',sp.simplify(u_s))
assert sp.simplify(u_s-3*v0*sp.sin(th)/4)==0
# momentum: horizontal 2mV = m*0 + m*V2x ; vertical 0 = -m u + m V2z
V2x = 2*V; V2z = u_s; print('(2)',V2x,V2z)
tB = [s for s in sp.solve(sp.Eq(h + V2z*t - g*t**2/2,0),t) if sp.simplify(s).subs({v0:1,th:0.7,g:1})>0][0]
tB=sp.simplify(tB); print('(3) tB=',tB); assert sp.simplify(tB-2*T)==0
xB2 = sp.simplify(V*T + V2x*tB); print('    x=',xB2); assert sp.simplify(xB2-5*V*T)==0
dK = sp.simplify(m*u_s**2/2 + m*(V2x**2+V2z**2)/2 - 2*m*V**2/2); print('(4) dK=',sp.expand(dK))
assert sp.simplify(dK - m*v0**2*(sp.cos(th)**2+sp.Rational(9,16)*sp.sin(th)**2))==0
# CM at t_A = T/2 lies on original parabola
xA1 = V*T; zA1=0; xA2 = V*T+V2x*T/2; zA2 = h+V2z*T/2-g*(T/2)**2/2
assert sp.simplify((xA1+xA2)/2-(V*T+V*T/2))==0 and sp.simplify((zA1+zA2)/2-(h-g*(T/2)**2/2))==0
xcm_end = (V*T + xB2)/2; d = sp.simplify(xcm_end - 2*V*T); print('(5) xcm-xB=',d)
# numeric: simulate both fragments from launch with explicit kinematics
v0v,thv,gv=20.0,np.radians(50),9.8
Vv=v0v*np.cos(thv); Tv=v0v*np.sin(thv)/gv; hv=gv*Tv**2/2; uv=0.75*v0v*np.sin(thv)
from scipy.optimize import brentq
t1=brentq(lambda s: hv-uv*s-gv*s*s/2,1e-9,100); t2=brentq(lambda s: hv+uv*s-gv*s*s/2,1e-9,100)
print('num t1/T',t1/Tv,' t2/T',t2/Tv,' x2/(VT)',(Vv*Tv+2*Vv*t2)/(Vv*Tv))
assert abs(t1/Tv-0.5)<1e-9 and abs(t2/Tv-2)<1e-9
print('ALL OK')
