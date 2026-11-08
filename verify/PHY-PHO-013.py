# 2000[1] 類題：光学的糖蜜（2本の対向レーザーによる摩擦力と反跳加熱のつり合い）
import sympy as sp
import numpy as np
h,c,nu0,M,v,R0,beta,kB,gam,dl=sp.symbols('h c nu0 M v R0 beta k_B gamma delta',positive=True)
p=h*nu0/c
# (1) 見かけの振動数（1次）
nu1=nu0*(1-v/c)   # +x向きの光（原子と同じ向き）
nu2=nu0*(1+v/c)   # -x向きの光（向かってくる）
R=lambda nup: R0*(1+beta*(nup-nu0))
F=sp.expand(p*R(nu1)-p*R(nu2))
alpha=sp.simplify(-F/v)
print('alpha =',alpha)
assert sp.simplify(alpha-2*h*nu0**2*R0*beta/c**2)==0
# 加熱：1回の散乱で2回のけり（吸収と放出），1回あたり <Δ(Mv)^2>=p^2，散乱率 2R0
heat=4*R0*p**2       # d<(Mv)^2>/dt の加熱分
v2=sp.symbols('v2',positive=True)
sol=sp.solve(sp.Eq(heat,2*M*alpha*v2),v2)[0]
K=sp.simplify(M*sol/2)
print('<K> =',K)
assert sp.simplify(K-h/(2*beta))==0
T=sp.simplify(2*K/kB)
print('T =',T)
# ローレンツ形の β
Delta=sp.symbols('Delta',real=True)
Rl=1/(1+Delta**2/gam**2)
b=sp.simplify((sp.diff(Rl,Delta)/Rl).subs(Delta,-dl))
print('beta(delta)=',b)
Td=sp.simplify(h/(kB*b))
dopt=sp.solve(sp.diff(Td,dl),dl)
print('delta_opt',dopt, 'Tmin',sp.simplify(Td.subs(dl,dopt[0])))
print('T numeric', 6.626e-34*4.9e6/1.381e-23)
# モンテカルロ（無次元）：p=1, M=1 として，各ステップで各ビームの吸収確率 ∝ R
rng=np.random.default_rng(1)
pp=1.0; Mm=50.0; r0=0.05; bb=0.0008   # beta*nu0*v/c の代わりに：R1=r0*(1-kappa*v), R2=r0*(1+kappa*v)
# 本文の R=R0[1+β(ν'-ν0)] で ν'-ν0=∓ν0 v/c。無次元化：kappa=β ν0/c，p=hν0/c → h/β = p/(kappa)·... 確認用に直接 κ を使う
kappa=0.02
N=4000; steps=150000
vel=np.zeros(N)
acc=[]
for s in range(steps):
    u=rng.random(N)
    a1=r0*(1-kappa*vel); a2=r0*(1+kappa*vel)
    ev1=u<a1; ev2=(u>=a1)&(u<a1+a2)
    ev=ev1|ev2
    kick=np.where(ev1,pp,0.0)+np.where(ev2,-pp,0.0)
    em=np.where(rng.random(N)<0.5,pp,-pp)*ev
    vel+= (kick+em)/Mm
    if s>60000 and s%100==0: acc.append(np.mean(vel**2))
Ksim=0.5*Mm*np.mean(acc)
# 理論：<K>=R0 p^2/alpha, alpha = 2 p kappa r0 (力 = p r0(-2 kappa v))
Kth=r0*pp**2/(2*pp*kappa*r0)
print('MC <K> =',Ksim,' theory =',Kth)
assert abs(Ksim/Kth-1)<0.05
