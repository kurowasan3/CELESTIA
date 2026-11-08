# P90-1: レーザー減速とゼーマン減速器（共鳴条件を保つ磁場分布）とエネルギー収支
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
h=6.63e-34; M=3.82e-26; lam=589e-9; c=3.00e8; V=6.0e2; nmax=3.0e7; kap=9.27e-24
dv=h/(M*lam); N=V/dv; amax=nmax*dv; a=amax/2; L=V**2/(2*a); B0=h*V/(kap*lam)
print('dv',dv,'N',N,'amax',amax,'L',L,'B0',B0)
# 記号: 共鳴 νL(1+v/c)=ν0+κB/h, νL=ν0 → B=h ν0 v/(κ c)=h v/(κ λ)
x,Vs,As=sp.symbols('x V a',positive=True)
v=sp.sqrt(Vs**2-2*As*x); print('B(x)/B0 =',sp.simplify(v/Vs))
# エネルギー収支: 吸収 hνL, 放出平均 hν0+κB, 差 = hνL v/c, 運動エネルギー減少 M v dv = h v/λ
nu0=sp.symbols('nu0',positive=True); vv,B,kk,hh,cc=sp.symbols('v B kappa h c',positive=True)
Bres=sp.solve(sp.Eq(nu0*(1+vv/cc),nu0+kk*B/hh),B)[0]
print('emit-absorb - h v/λ =',sp.simplify((hh*nu0+kk*Bres-hh*nu0)-hh*vv*nu0/cc))
# シミュレーション: 原子は |共鳴ずれ| < 許容幅 のときのみ率 nmax で減速（離散キックを連続化）
tol=0.5  # m/s 相当の速度許容幅
def run(vin):
    def f(t,y):
        xx,vv_=y
        vres=V*np.sqrt(max(0,1-xx/L)) if xx<L else 0
        acc=-amax if abs(vv_-vres)<tol else 0
        return [vv_,acc]
    s=solve_ivp(f,[0,0.05],[0,vin],max_step=2e-6,events=lambda t,y:y[0]-1.2*L)
    return s
for vin in [600.0,450.0,650.0]:
    s=run(vin); xs,vs=s.y
    # 減速が始まる位置
    idx=np.argmax(vs<vin-1.0) if (vs<vin-1.0).any() else None
    print('vin',vin,'start decel x',xs[idx] if idx else None,'theory',L*(1-(vin/V)**2) if vin<=V else 'never','final v at L',np.interp(L*0.999,xs,vs))
