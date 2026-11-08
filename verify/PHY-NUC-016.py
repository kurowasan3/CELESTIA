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
# 数値
k0=9.0e9; e=1.6e-19; Z=79; RNv=7.0e-15
dK=lambda KeV: 2*k0*Z*e/(KeV)   # K in eV -> d [m]
print('d(5MeV)=',dK(5.0e6))
print('U(RN) MeV=',2*k0*Z*e/RNv/1e6)
print('K for sigma=half MeV=',4*k0*Z*e/RNv/1e6)
K40=40e6; dd=dK(K40); print('sigma ratio at 40MeV:',1-dd/RNv, ' b_max=',RNv*np.sqrt(1-dd/RNv))
