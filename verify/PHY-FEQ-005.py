# P10 類題：角度 θ で交わる2枚の鉛直な壁の隅で体を支える（摩擦の水平成分が必要）
import sympy as sp, numpy as np
from scipy.optimize import minimize
th,N1,N2,m,g,mu=sp.symbols('theta N1 N2 m g mu',positive=True)
a,b=sp.symbols('a b')
e1=sp.Matrix([1,0]); e2=sp.Matrix([sp.cos(th),sp.sin(th)])
n1=sp.Matrix([0,1]); n2=sp.Matrix([sp.sin(th),-sp.cos(th)])
sol=sp.solve(list(N1*n1+N2*n2+a*e1+b*e2),[a,b])
print('a=',sp.simplify(sol[a]),' b=',sp.simplify(sol[b]))
# 対称 N1=N2=N
N=sp.symbols('N',positive=True)
asym=sp.simplify(sol[a].subs({N1:N,N2:N})); print('symmetric horizontal friction =',sp.simplify(asym), ' -> ', sp.simplify(asym+N*sp.tan(th/2)))
Nmin=m*g/(2*sp.sqrt(mu**2-sp.tan(th/2)**2))
print('Nmin=',Nmin, float(Nmin.subs({m:60,g:9.8,mu:0.8,th:sp.pi/3})))
print('theta max for mu=0.8:',float(2*sp.atan(0.8)*180/sp.pi))
# 非対称でも最小必要μ = tan(θ/2) を数値で確認：与えた θ で、比 r=N2/N1 と鉛直分担を自由にして必要 μ を最小化
def need_mu(params,thv):
    r,s=params  # r=N2/N1, s = fraction of weight on wall1
    if r<=0 or not (0<=s<=1): return 1e9
    N1v=1.0; N2v=r
    av=float(sol[a].subs({N1:N1v,N2:N2v,th:thv})); bv=float(sol[b].subs({N1:N1v,N2:N2v,th:thv}))
    W=1e-9  # 重さ→0 の極限（N→∞）で必要μが最小
    mu1=np.hypot(av,s*W)/N1v; mu2=np.hypot(bv,(1-s)*W)/N2v
    return max(mu1,mu2)
for deg in [40,60,90,120]:
    thv=np.radians(deg)
    best=min((minimize(need_mu,[r0,0.5],args=(thv,),method='Nelder-Mead').fun for r0 in [0.3,0.7,1,1.5,3]))
    print(deg,'min mu numeric',best,' tan(th/2)=',np.tan(thv/2))
    assert abs(best-np.tan(thv/2))<1e-4
# 有限の重さでの最小 N（対称・非対称を探索）: θ=60°, μ=0.8, m=60
thv=np.pi/3; muv=0.8; Wt=60*9.8
def feasible(N1v,N2v,s):
    av=float(sol[a].subs({N1:N1v,N2:N2v,th:thv})); bv=float(sol[b].subs({N1:N1v,N2:N2v,th:thv}))
    return np.hypot(av,s*Wt)<=muv*N1v and np.hypot(bv,(1-s)*Wt)<=muv*N2v
# 「両壁を押す力の大きい方」を最小化（グリッド）
best=1e9
for N1v in np.linspace(400,700,151):
    for r in np.linspace(0.8,1.25,46):
        for s in np.linspace(0.3,0.7,41):
            if feasible(N1v,N1v*r,s): best=min(best,max(N1v,N1v*r))
print('grid min of max normal force',best,' theory',float(Nmin.subs({m:60,g:9.8,mu:0.8,th:sp.pi/3})))
