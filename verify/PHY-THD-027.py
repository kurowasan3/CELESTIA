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
# ---- 新しい (6)：一般の比熱比 γ で，吸熱 0 の高さと安定限界が一致する
gm=sp.symbols('gamma',positive=True)
dQg=n*R*sp.diff(T,y)/(gm-1)+P*S
yq=sp.solve(sp.factor(dQg),y); print('(6) 吸熱0:',yq)
dFg=rho*g*S-gm*P*S/y
ys=sp.solve(dFg,y); print('    安定限界:',ys)
ycg=gm*(a+H)/(gm+1)
assert len(yq)==1 and len(ys)==1 and sp.simplify(yq[0]-ycg)==0 and sp.simplify(ys[0]-ycg)==0
assert sp.simplify(ycg.subs(gm,sp.Rational(5,3))-yc)==0
width=sp.simplify(ycg-ym); print('    幅',sp.factor(width),' γ→1:',sp.limit(width,gm,1))
assert sp.simplify(width-(gm-1)*(a+H)/(2*(gm+1)))==0 and sp.limit(width,gm,1)==0
# 数値：γ=7/5, a=0.2, H=1 で y_c の前後でモル比熱の符号と断熱変位への応答を調べる
gv=1.4; an,Hn=0.2,1.0; ycv=gv*(an+Hn)/(gv+1)
Cf=sp.lambdify(y,sp.simplify((dQg/(n*sp.diff(T,y))).subs({gm:gv,a:an,H:Hn})))
print('    γ=1.4: yc=',ycv,' C/R at ym+0.01:',Cf(0.61)/1,' at yc+0.01:',Cf(ycv+0.01))
for y0f,expect in ((ycv-0.02,'stable'),(ycv+0.02,'unstable')):
    P0g=an+Hn-y0f
    def rhs(t,s_):
        yy,v=s_; Pg=P0g*(y0f/yy)**gv
        liq=(an+Hn-yy) if yy>=y0f else (an+Hn-y0f)
        return [v,(Pg-liq)/1e-3]
    so=solve_ivp(rhs,[0,1],[y0f+1e-4,0],max_step=1e-3)
    md=np.max(np.abs(so.y[0]-y0f)); print('    y0=',round(y0f,3),expect,' max disp',md)
    assert (md<1e-3) if expect=='stable' else (md>1e-2)
print('OK (6)')
