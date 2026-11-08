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
# (2) の N_min を有限の重さで数値確認（θ=60°, μ=0.8, m=60 は検算用の値）
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

# ---- (4) different mu1, mu2: feasible (for large push) iff atan(mu1)+atan(mu2) > theta ----
# weightless limit: N1=1, N2=r; need |a|/N1 < mu1 and |b|/N2 < mu2 for some r>0
fa=sp.lambdify((N1,N2,th),sol[a]); fb=sp.lambdify((N1,N2,th),sol[b])
def margin(mu1,mu2,thv):
    rs=np.exp(np.linspace(-6,6,24001))
    m1=np.abs(fa(1.0,rs,thv))/1.0/mu1; m2=np.abs(fb(1.0,rs,thv))/rs/mu2
    return np.min(np.maximum(m1,m2))   # <1 feasible
rng=np.random.default_rng(0); bad=0; cnt=0
for _ in range(3000):
    mu1,mu2=rng.uniform(0.05,3,2); thv=rng.uniform(0.05,np.pi-0.05)
    crit=np.arctan(mu1)+np.arctan(mu2)-thv
    if abs(crit)<5e-3: continue
    cnt+=1
    if (margin(mu1,mu2,thv)<1)!=(crit>0): bad+=1
print('(4) random tests',cnt,'mismatches',bad); assert bad==0
# theta=90deg: atan(mu1)+atan(mu2)>pi/2 <=> mu1*mu2>1
for mu1,mu2 in [(0.5,2.1),(0.5,1.9),(1.2,0.85),(3,0.3)]:
    assert (np.arctan(mu1)+np.arctan(mu2)>np.pi/2)==(mu1*mu2>1)
print('(4) right angle: mu1*mu2>1 equivalence OK')
# sufficiency with finite weight: construct K large, phi split, check friction inequality
mu1,mu2,thv,W=0.8,1.5,np.radians(90),1.0
l1,l2=np.arctan(mu1),np.arctan(mu2); assert l1+l2>thv
ph1=thv*l1/(l1+l2); ph2=thv-ph1
K=20.0
ok=all(np.hypot(K*np.sin(p_),W/2)<=mu_*K*np.cos(p_) for p_,mu_ in [(ph1,mu1),(ph2,mu2)])
print('(4) construction with K=20 feasible:',ok); assert ok
