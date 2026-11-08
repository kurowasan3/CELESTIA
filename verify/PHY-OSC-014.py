# M42-1: column (cross-section S, height h, density 2rho/3) in container cross-section 3S.
# released with bottom touching water surface. y = descent of column (absolute).
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
rho,S,h,g,y=sp.symbols('rho S h g y',positive=True)
m=sp.Rational(2,3)*rho*S*h
# water volume conservation: submerged depth d, level rise H: 3S*H = S*d, d = y + H
H,d=sp.symbols('H d')
sol=sp.solve([sp.Eq(3*S*H,S*d),sp.Eq(d,y+H)],[H,d]); print('level rise, depth:',sol)
dd=sol[d]
F=m*g-rho*S*dd*g
ye=sp.solve(sp.Eq(F,0),y)[0]; print('(1) y_eq',ye,' depth',dd.subs(y,ye))
K=-sp.diff(F,y); w2=sp.simplify(K/m); print('(2) omega^2',w2,' T=',sp.simplify(2*sp.pi/sp.sqrt(w2)))
yf=sp.solve(sp.Eq(dd,h),y)[0]; print('top flush at y=',yf)
A=ye; v2=w2*(A**2-(yf-ye)**2); print('(3) v^2=',sp.simplify(v2))
a=sp.simplify((m*g-rho*S*h*g)/m); print('   accel after full submersion',a)
dmax=sp.simplify(v2/(2*(-a))); print('(4) extra depth',dmax)
tS=sp.acos(-(yf-ye)/A)/sp.sqrt(w2)   # y-ye=-A cos(wt)
print('   t_SHM',sp.simplify(tS))
td=sp.sqrt(v2)/(-a); P=2*(tS+td); print('(5) period',sp.simplify(P))
# numeric
rr,SS,hh,gg=1000.0,0.01,0.3,9.8
mm=2/3*rr*SS*hh
def f(t,Y):
    yy,vv=Y; dep=min(1.5*yy,hh) if yy>0 else 0.0
    return [vv,(mm*gg-rr*SS*dep*gg)/mm]
ev=lambda t,Y: Y[1]; ev.direction=-1; ev.terminal=False
r=solve_ivp(f,[0,3],[0,0],rtol=1e-11,atol=1e-13,max_step=1e-4,events=ev)
sub={rho:rr,S:SS,h:hh,g:gg}
print('num ymax',r.y[0].max(),' th', float((yf+dmax).subs(sub)))
print('num period (v=0 return times)', r.t_events[0][:3], ' th', float(P.subs(sub)))
# depth of top below surface at max: column descends extra dmax with level fixed
