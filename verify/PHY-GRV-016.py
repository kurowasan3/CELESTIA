# Q087-1: 正三角形・正方形の頂点の質点の重力による崩壊と回転（有効質量と細長い楕円）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
G,m,R=sp.symbols('G m R',positive=True)
def Meff(n):
    # 中心向きの力 = G m * Meff / R^2
    tot=0
    for j in range(1,n):
        ang=2*sp.pi*j/n; dist=2*R*sp.sin(ang/2)
        tot+=G*m*m/dist**2*sp.sin(ang/2)   # 中心向き成分 = cos(半角の余角) = sin(ang/2)
    return sp.nsimplify(sp.simplify(tot*R**2/(G*m)))
M3=sp.simplify(Meff(3)); M4=sp.simplify(Meff(4))
print('Meff3=',M3,' Meff4=',M4)
assert sp.simplify(M3-m/sp.sqrt(3))==0 and sp.simplify(M4-m*(1+2*sp.sqrt(2))/4)==0
tf=lambda Me: sp.pi*sp.sqrt(R**3/(8*G*Me))
t3=tf(M3); t4=tf(M4); print('t3=',sp.simplify(t3),' t4/t3=',sp.nsimplify(sp.simplify(t4/t3)),float(t4/t3))
# 数値：n体を直接積分
def nbody(n,Gv=1.0,mv=1.0,Rv=1.0,vt=0.0,tend=10.0,stop=0.01):
    ang=2*np.pi*np.arange(n)/n
    pos=np.stack([Rv*np.cos(ang),Rv*np.sin(ang)],1); vel=vt*np.stack([-np.sin(ang),np.cos(ang)],1)
    def f(t,z):
        p=z[:2*n].reshape(n,2); v=z[2*n:].reshape(n,2); a=np.zeros_like(p)
        for i in range(n):
            for j in range(n):
                if i!=j:
                    dr=p[j]-p[i]; a[i]+=Gv*mv*dr/np.linalg.norm(dr)**3
        return np.concatenate([v.ravel(),a.ravel()])
    ev=lambda t,z: np.linalg.norm(z[:2])-stop; ev.terminal=True
    return solve_ivp(f,[0,tend],np.concatenate([pos.ravel(),vel.ravel()]),events=ev,rtol=1e-11,atol=1e-13,dense_output=True)
for n,Me in [(3,M3),(4,M4)]:
    s=nbody(n,stop=1e-4); th=float(tf(Me).subs({G:1,m:1,R:1}))
    print('n=',n,' t_num(r->1e-4)=',s.t[-1],' t_th=',th); assert abs(s.t[-1]-th)<1e-3
# 回転：v0^2 = G Meff/(2R) のとき最小距離 R/3、そこまでの時間
vc2=G*M3/R; v0=sp.sqrt(vc2/2)
a=sp.Symbol('a')
rmin=R*(v0**2*R/(G*M3))/(2-v0**2*R/(G*M3)); print('rmin=',sp.simplify(rmin))
half=sp.pi*sp.sqrt(((R+rmin)/2)**3/(G*M3)); print('time to rmin=',sp.simplify(half))
s=nbody(3,vt=float(v0.subs({G:1,m:1,R:1})),tend=float(half.subs({G:1,m:1,R:1}))*1.5,stop=-1)
ts=np.linspace(0,s.t[-1],300001); r=np.hypot(s.sol(ts)[0],s.sol(ts)[1]); i=np.argmin(r)
print('num rmin=',r[i],' at t=',ts[i],' th',1/3,float(half.subs({G:1,m:1,R:1})))
assert abs(r[i]-1/3)<1e-5 and abs(ts[i]-float(half.subs({G:1,m:1,R:1})))<1e-3
print('min side =',sp.sqrt(3)*sp.Rational(1,3),'R')
print('OK')
