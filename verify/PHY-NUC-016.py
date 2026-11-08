# J148-1 検証: 衝突径数をもつα粒子の最接近距離と原子核に触れる断面積
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
r,b,d,K=sp.symbols('r b d K',positive=True)
# エネルギー: K = K b^2/r^2 + K d/r
rmin=sp.solve(sp.Eq(K, K*b**2/r**2 + K*d/r), r); print('rmin',rmin)
rm=[s for s in rmin if s.subs({b:1,d:1})>0][0]
print('b for rmin=2d:',sp.solve(sp.Eq(rm,2*d),b))
RN=sp.symbols('R_N',positive=True)
print('b_max^2:',sp.solve(sp.Eq(rm,RN),b))
# 数値軌道: m=1, K=1(v0=sqrt2), d=1 -> k*Q*q = K*d = 1
def num_rmin(bb):
    m=1;C=1.0;r0=np.hypot(20000,bb);v0=np.sqrt(2*(1-C/r0))  # 全エネルギーをちょうど K=1 にそろえる
    def f(t,y):
        rr=np.hypot(y[0],y[1]); return [y[2],y[3],C*y[0]/rr**3,C*y[1]/rr**3]
    y0=[-20000,bb,v0,0]
    s=solve_ivp(f,[0,30000],y0,rtol=1e-10,atol=1e-12,max_step=1.0, dense_output=True)
    tt=np.linspace(0,30000,3000001); xy=s.sol(tt); return np.min(np.hypot(xy[0],xy[1]))
    return np.min(np.hypot(s.y[0],s.y[1]))
for bb in [0.0,0.5,np.sqrt(2),3.0]:
    print('b',bb,'num',num_rmin(bb),'theory',(1+np.sqrt(1+4*bb**2))/2)

# ---- new (5): finite nuclear mass (recoil) ----
m_,M_,v0,UN=sp.symbols('m M v_0 U_N',positive=True)
V=m_*v0/(m_+M_); Kv=m_*v0**2/2
print('KE at closest (head-on) / K =',sp.simplify((m_+M_)*V**2/2/Kv))
print('relative KE 1/2 mu v0^2 / K =',sp.simplify(m_*M_/(m_+M_)*v0**2/2/Kv))
sigma=sp.pi*RN**2*(1-(1+m_/M_)*UN/K)
print('sigma M->oo:',sp.limit(sigma,M_,sp.oo))
# numerical two-body check: nucleus free (mass M), alpha mass m; find minimum separation vs theory
def two_body_rmin(mm,MM,bb,C=1.0,Kin=1.0):
    v0=np.sqrt(2*Kin/mm); X=4000.0
    def f(t,y):
        dx=y[0]-y[2]; dy=y[1]-y[3]; rr=np.hypot(dx,dy); F=C/rr**3
        return [y[4],y[5],y[6],y[7], F*dx/mm, F*dy/mm, -F*dx/MM, -F*dy/MM]
    y0=[-X,bb,0,0,v0,0,0,0]
    s=solve_ivp(f,[0,2*X/v0],y0,rtol=1e-10,atol=1e-12,max_step=0.5,dense_output=True)
    tt=np.linspace(0,2*X/v0,1000001); Y=s.sol(tt); return np.min(np.hypot(Y[0]-Y[2],Y[1]-Y[3]))
mm,MM=1.0,3.0; Kin=1.0; C=1.0
Krel=MM/(mm+MM)*Kin; dp=C/Krel
for bb in [0.0,1.0,2.5]:
    print(' two-body b=',bb,' num rmin',two_body_rmin(mm,MM,bb),' theory',(dp+np.sqrt(dp**2+4*bb**2))/2)
