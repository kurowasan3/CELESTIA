# P30-1 平均自由行程と熱伝導（連続領域とクヌーセン領域）
import sympy as sp
k,T,p,d,vb,L,dT=sp.symbols('k T p d vbar L DeltaT',positive=True)
nden=p/(k*T)
lam=sp.solve(sp.Eq(nden*sp.pi*d**2*sp.Symbol('lam'),1),sp.Symbol('lam'))[0]; print('lambda=',lam)
f=sp.Rational(1,4)*nden*vb
# model A: molecules come from distance lam on each side, carry 3/2 k T
grad=dT/L
qA=sp.simplify(f*sp.Rational(3,2)*k*(2*lam*grad)); print('qA=',qA,' p-dependence:',sp.diff(qA,p))
qB=sp.simplify(f*sp.Rational(3,2)*k*dT); print('qB=',qB)
ratio=sp.simplify(qB/qA); print('qB/qA=',sp.simplify(ratio), ' check L/(2lam):',sp.simplify(ratio-L/(2*lam)))
pc=sp.solve(sp.Eq(qA,qB),p)[0]; print('crossover p=',pc)
# 参考（数値例，問題文には含まない）: lambda0=6.8e-8 at 1.0e5 Pa; L=5.0e-3; want q = qA/100 in B regime -> L/(2 lam)=1/100
lam0,p0,Lv=6.8e-8,1.0e5,5.0e-3
lam_needed=50*Lv; pneed=p0*lam0/lam_needed; print('p needed=',pneed)
# lambda at crossover numeric
print('p_c numeric =',p0*lam0/(Lv/2))
# check 1/100 by direct formulas with Ar-like numbers
import math
kB=1.38e-23;Tv=300; dv=math.sqrt(kB*Tv/(math.pi*p0*lam0)); vbv=400
qAv=float(qA.subs({k:kB,T:Tv,d:dv,vb:vbv,L:Lv,dT:1}))
qBv=float(qB.subs({k:kB,T:Tv,p:pneed,vb:vbv,dT:1}))
print('d=',dv,'qB/qA=',qBv/qAv)
# flux relative to 1 atm at 1.0e2 Pa and 1.0 Pa (actual flux = min(qA,qB))
for pv in [1.0e2,1.0]:
    lamv=lam0*p0/pv
    rel=min(1.0, Lv/(2*lamv))
    print('p=',pv,'lambda=',lamv,'q/q_atm=',rel)
# ---- 新しい (5)：N 枚の板を入れたときの定常熱流（モデルB／モデルA）
for Nn in range(1,6):
    taus=sp.symbols('tau1:%d'%(Nn+1)); tt=[T+dT]+list(taus)+[T]
    for model in ('B','A'):
        gap=L/(Nn+1)
        flux=[ (sp.Rational(3,8)*p*vb/T*(tt[i]-tt[i+1])) if model=='B'
               else (3*k*vb/(4*sp.pi*d**2)*(tt[i]-tt[i+1])/gap) for i in range(Nn+1)]
        sol=sp.solve([sp.Eq(flux[i],flux[i+1]) for i in range(Nn)],taus,dict=True)[0]
        q=sp.simplify(flux[0].subs(sol))
        ref=qB if model=='B' else qA
        r_=sp.simplify(q/ref)
        assert r_==(sp.Rational(1,Nn+1) if model=='B' else 1),(Nn,model,r_)
        # 入れた板の温度は等間隔に下がる
        assert all(sp.simplify(sol[taus[j]]-(T+dT*(Nn-j)/(Nn+1)))==0 for j in range(Nn))
print('(5) モデルB: 1/(N+1) 倍，モデルA: 1 倍（N=1..5）を確認')
# 粒子的な独立検証（モデルB）：板の間を飛ぶ分子が板の温度をもってはね返る1次元の鎖を，
# 各板のエネルギー収支 dE/dt = 0 となるまで時間発展させ，定常の熱流を比べる
import numpy as np
def chain(Nn,steps=200000,dt=0.01):
    Tm=np.linspace(1.0,0.0,Nn+2)[1:-1].copy()*0+0.5   # 初期温度は適当
    for _ in range(steps):
        full=np.concatenate(([1.0],Tm,[0.0]))
        fl=full[:-1]-full[1:]          # 各すき間の流れ（係数1）
        Tm+=dt*(fl[:-1]-fl[1:])
    full=np.concatenate(([1.0],Tm,[0.0])); return (full[:-1]-full[1:]).mean()
for Nn in (1,3): print('chain N=',Nn,' flux/flux0=',chain(Nn)); assert abs(chain(Nn)-1/(Nn+1))<1e-6
print('OK (5)')
