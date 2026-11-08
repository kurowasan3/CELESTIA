# O20211-1 : free tether co-rotating with the Earth; tension maximum; release from the tether
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
r,R0,Rs,lam,w,G,M,rp,rr=sp.symbols('r R_0 R_s lambda omega G M r_p rr',positive=True)
GM=w**2*Rs**3
# element: lam dr * w^2 r = lam dr GM/r^2 - dT  -> dT/dr = lam(GM/r^2 - w^2 r), T(R0)=0
T=sp.integrate(lam*(GM/rr**2-w**2*rr),(rr,R0,r))
print('T(r)=',sp.simplify(T))
print('dT/dr=0 at',sp.solve(sp.diff(T,r),r))
Tmax=sp.factor(sp.simplify(T.subs(r,Rs))); print('Tmax=',Tmax)
# top end: T(R2)=0
R2=sp.symbols('R_2',positive=True)
print('R2 eq:',sp.factor(sp.simplify(T.subs(r,R2))))
# release at r with speed w r (tangential). other apsis from L & E conservation
v=w*r
sol=sp.solve([sp.Eq(r*v,rp*sp.Symbol('vp')), sp.Eq(v**2/2-GM/r, sp.Symbol('vp')**2/2-GM/rp)],[rp,sp.Symbol('vp')],dict=True)
for s_ in sol: print('other apsis:',sp.simplify(s_[rp]))
s=sp.symbols('s',positive=True)
print('in s=(r/Rs)^3: r*s/(2-s) check:',sp.simplify(r*(r/Rs)**3/(2-(r/Rs)**3)-[q[rp] for q in sol if sp.simplify(q[rp]-r)!=0][0]))
# orbit integration check of (4) (units R0=1, GM=1, Rs=6.6)
GMv=1.0; Rsv=6.6; wv=np.sqrt(GMv/Rsv**3)
def peri(r0):
    def f(t,y):
        x,yy,vx,vy=y; d=np.hypot(x,yy)**3; return [vx,vy,-GMv*x/d,-GMv*yy/d]
    s=solve_ivp(f,[0,2e3],[r0,0,0,wv*r0],rtol=1e-11,atol=1e-12,max_step=0.5)
    return np.min(np.hypot(s.y[0],s.y[1])), np.max(np.hypot(s.y[0],s.y[1]))
for r0 in [4.0,5.5,7.5]:
    s3=(r0/Rsv)**3; print('r0=%.3f numeric min/max r'%r0,peri(r0),' formula other apsis',r0*s3/(2-s3))
# ---- (5) top end, total gravity, centre of mass ----
R2s=(-R0+sp.sqrt(R0**2+8*Rs**3/R0))/2
assert sp.simplify(sp.expand(T.subs(r,R2s)))==0
print('(5) R2 =',R2s,' satisfies T(R2)=0')
F=sp.integrate(G*M*lam/rr**2,(rr,R0,R2)); mt=lam*(R2-R0)
assert sp.simplify(F-G*M*mt/(R0*R2))==0
rG=(R0+R2)/2
# whole-tether equation m rG w^2 = F equivalent to T(R2)=0
print('(5) m rG w^2 - F with GM=w^2Rs^3:',sp.factor(sp.simplify(mt*rG*w**2-lam*GM*(1/R0-1/R2))))
print('    T(R2)*?:',sp.factor(sp.simplify(T.subs(r,R2))))
for ratio in [1.5,2.0,6.6,20.0]:
    R0v=1.0; Rsv_=ratio; R2v=float(R2s.subs({R0:R0v,Rs:Rsv_}))
    rGv=(R0v+R2v)/2
    print('    Rs/R0=%.1f R2=%.4f rG=%.4f rG>Rs:%s  Rs^3/(R0R2)=%.4f'%(ratio,R2v,rGv,rGv>Rsv_,Rsv_**3/(R0v*R2v)))
    assert rGv>Rsv_ and abs(rGv-Rsv_**3/(R0v*R2v))<1e-9
# numeric force sum
R0v,Rsv_=1.0,6.6; R2v=float(R2s.subs({R0:R0v,Rs:Rsv_})); xs=np.linspace(R0v,R2v,200001)
Fg=np.trapezoid(1/xs**2,xs)*Rsv_**3   # GM=Rs^3 (w=1), lambda=1
Fc=np.trapezoid(xs,xs)
print('(5) numeric total gravity %.6f vs total centripetal %.6f'%(Fg,Fc)); assert abs(Fg-Fc)<1e-6*Fc
