# K20046-1 verification: stiff elastic string = elastic circular wall; retracing projectile orbit
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,l,y0,u0,T=sp.symbols('g l y0 u0 T',real=True)
x1=u0*T; y1=y0-g*T**2/2; vx=u0; vy=-g*T
cond=sp.simplify(vy*x1-vx*y1)  # velocity parallel to OD
print('parallel condition:',sp.factor(cond),' -> T^2 =',sp.solve(cond,T**2) if False else sp.solve(sp.Eq(cond,0),T))
T2=-2*y0/g
eq=sp.Eq((u0**2*T2)+(y0-g*T2/2)**2,l**2)
print('u0^2 =',sp.solve(eq,u0**2))
u0sq=sp.solve(eq,u0**2)[0]
print('y0=-l/4: u0^2=',sp.simplify(u0sq.subs(y0,-l/4)),' D=',sp.simplify((sp.sqrt(u0sq*T2)).subs(y0,-l/4)),sp.simplify((y0-g*T2/2).subs(y0,-l/4)))
# numeric with finite k
gn,ln,mn=9.8,1.0,1.0
y0n=-ln/4; u0n=float(sp.sqrt(u0sq).subs({g:gn,l:ln,y0:y0n})); Tn=np.sqrt(-2*y0n/gn)
for kn in (1e5,1e7):
    def rhs(t,s):
        x,y,vx,vy=s; r=np.hypot(x,y); fx=fy=0.0
        if r>ln:
            f=-kn*(r-ln); fx=f*x/r; fy=f*y/r
        return [vx,vy,fx/mn,fy/mn-gn]
    sol=solve_ivp(rhs,[0,2*Tn+0.05],[0,y0n,u0n,0],rtol=1e-11,atol=1e-12,max_step=1e-4,dense_output=True)
    # time of return to x=0 after bounce
    ts=np.linspace(Tn,2*Tn+0.05,200001); xs=sol.sol(ts)[0]
    j=np.where((xs[:-1]>0)&(xs[1:]<=0))[0][0]; tr=ts[j]; st=sol.sol(tr)
    smax=max(np.hypot(*sol.sol(t)[:2])-ln for t in np.linspace(Tn-0.01,Tn+0.02,30001))
    print(f'k={kn:.0e}: return to x=0 at t={tr:.5f} (2T={2*Tn:.5f}) y={st[1]:.5f} (y0={y0n}) v=({st[2]:.4f},{st[3]:.4f}) vs (-{u0n:.4f},0); max stretch={smax:.5f} est sqrt(2glm/k)={np.sqrt(2*gn*ln*mn/kn):.5f}')
print('period 4T =',4*Tn,' = 2 sqrt(2l/g)=',2*np.sqrt(2*ln/gn))
# existence range: u0^2>0 requires l^2-4y0^2>0 with y0<0
print('range: ', sp.solve(sp.And(u0sq.subs({g:1,l:1})>0), y0) if False else 'u0^2 = g(l^2-4y0^2)/(-2y0) >0 -> -l/2<y0<0')
