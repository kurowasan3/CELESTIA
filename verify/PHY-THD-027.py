import sympy as sp
rho,g,S,H,a,n,R,y=sp.symbols('rho g S H a n R y',positive=True)
P=rho*g*(a+H-y); T=P*S*y/(n*R)
ym=sp.solve(sp.diff(T,y),y)[0]; print('ym',ym,' Tmax',sp.factor(T.subs(y,ym)))
dQdy=sp.Rational(3,2)*n*R*sp.diff(T,y)+P*S
print('dQ/dy',sp.factor(dQdy))
C=sp.simplify(dQdy/(n*sp.diff(T,y))); print('C',sp.factor(C))
yc=sp.solve(dQdy,y)[0]; print('yc',yc)
# adiabatic stiffness: force up on piston F=P_gas S - (P0+rho g (H-y))S ; dP_gas/dy adiabatic = -5/3 P/y
dF=(-sp.Rational(5,3)*P/y)*S + rho*g*S
print('stability boundary',sp.solve(dF,y))
sub={a:H/5}
Q=sp.integrate(dQdy.subs(sub),(y,H/5,sp.Rational(3,4)*H)); print('Q',sp.factor(Q))
print('T(yc)/Tmax',sp.simplify((T.subs(y,yc)/T.subs(y,ym)).subs(sub)))
print('ym,yc with a=H/5',ym.subs(sub),yc.subs(sub))
# first law check: Q = dU + W
W=sp.integrate((P*S).subs(sub),(y,H/5,sp.Rational(3,4)*H)); dU=sp.Rational(3,2)*n*R*(T.subs(sub).subs(y,sp.Rational(3,4)*H)-T.subs(sub).subs(y,H/5))
print('W',sp.factor(W),' dU',sp.factor(dU),' sum',sp.factor(W+dU))
# numeric stability check: small displacement adiabatic simulation at y slightly below/above yc
import numpy as np
from scipy.integrate import solve_ivp
for y0f in (0.70,0.80):
    Hn=1.0; an=0.2; rg=1.0; Sn=1.0; mp=1e-3
    y0=y0f; P0g=rg*(an+Hn-y0)
    def rhs(t,s):
        yy,v=s; Pg=P0g*(y0/yy)**(5/3); liq=rg*(an+Hn-yy)  # liquid always to rim (overflow when going up; for down assume refilled? only test up)
        return [v,(Pg-liq)*Sn/mp]
    sol=solve_ivp(rhs,[0,0.05],[y0+1e-4,0],max_step=1e-5)
    print('y0',y0f,'disp after',sol.y[0][-1]-y0, ' (yc=0.75)')
for y0f in (0.70,0.74,0.76,0.80):
    y0=y0f; P0g=1.0*(1.2-y0)
    def rhs(t,s):
        yy,v=s; Pg=P0g*(y0/yy)**(5/3); top=max(yy,y0)  # liquid surface: rim if piston above y0 (overflowed), else stays
        liq=1.0*(0.2+1.0-yy) if yy>=y0 else 1.0*(0.2+1.0-y0)
        return [v,(Pg-liq)/1e-3]
    sol=solve_ivp(rhs,[0,2],[y0+1e-4,0],max_step=1e-3)
    print('y0',y0f,'max disp',np.max(np.abs(sol.y[0]-y0)))
