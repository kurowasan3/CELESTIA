# P28-1 浮沈子のヒステリシス検証
import sympy as sp
n,R,T,m,rho,g,P,P0,z,H,d=sp.symbols('n R T m rho g P P0 z H delta',positive=True)
V=n*R*T/(P+rho*g*z)
Fup=rho*g*V-m*g
# (1) lid z=0 : leave lid when Fup<0
P1=sp.solve(sp.Eq(Fup.subs(z,0),0),P)[0]; print('P1=',P1)
# (2) equilibrium depth at P0
zc=sp.solve(sp.Eq(Fup.subs(P,P0),0),z)[0]; print('zc=',sp.simplify(zc))
# (3) net downward force at zc+delta (linear)
F=-Fup.subs(P,P0)
k=sp.simplify(sp.diff(F,z).subs(z,zc)); print('dF_down/dz=',k)
# (4) P2 at bottom
P2=sp.solve(sp.Eq(Fup.subs(z,H),0),P)[0]; print('P2=',sp.simplify(P2),' P1-P2=',sp.simplify(P1-P2))
Hmax=sp.solve(sp.Eq(P2,P0/2),H)[0]; print('Hmax=',sp.simplify(Hmax))
# numeric sanity: simulate sinking from lid with P slightly above P1 then reset to P0 at depth z0
import numpy as np
from scipy.integrate import solve_ivp
vals={n:1e-4,R:8.31,T:300,m:2.0e-3,rho:1000,g:9.8,P0:1.0e5}
Pv1=float(P1.subs(vals)); zcv=float(zc.subs(vals)); print('num P1',Pv1,'zc',zcv)
def run(z0,Pval):
    Ff=sp.lambdify((z,),Fup.subs(vals).subs(P,Pval))
    f=lambda t,y:[y[1], -float(Ff(y[0]))/vals[m] - 0.5*y[1]]
    top=lambda t,y:y[0]; top.terminal=True
    bot=lambda t,y:y[0]-10; bot.terminal=True
    s=solve_ivp(f,[0,600],[z0,0],max_step=0.05,events=[top,bot]); return round(s.y[0,-1],3)
print('reset at zc-0.05 ->',run(zcv-0.05,vals[P0]),' reset at zc+0.05 ->',run(zcv+0.05,vals[P0]))
# overdamped motion (inertia negligible): dz/dt = -Fup/gamma, stop at lid/bottom
def run_od(z0,Pval):
    Ff=sp.lambdify((z,),Fup.subs(vals).subs(P,Pval))
    top=lambda t,y:y[0]; top.terminal=True
    bot=lambda t,y:y[0]-10; bot.terminal=True
    s=solve_ivp(lambda t,y:[-float(Ff(y[0]))/0.01],[0,1e4],[z0],events=[top,bot],max_step=1)
    return round(s.y[0,-1],3)
print('overdamped: reset at zc-0.01 ->',run_od(zcv-0.01,vals[P0]),' zc+0.01 ->',run_od(zcv+0.01,vals[P0]))
P2v=float(P2.subs(vals).subs(H,10)); print('P2(H=10)=',P2v,' just below ->',run_od(10-1e-6,P2v-10),' just above ->',run_od(10-1e-6,P2v+10))
