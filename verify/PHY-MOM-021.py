# 2002[1] 類題：球を1個ずつ後方へ打ち出す台車（一度に打ち出す場合・連続噴射の極限との比較）
import sympy as sp, numpy as np
from fractions import Fraction as Fr
m,M,w,W,Wp=sp.symbols('m M w W Wp',positive=True)
# 1回の打ち出し（打ち出し後の相対速さ w）：質量 Mk-1 → Mk = Mk-1 - m
Mk1=sp.symbols('Mk1',positive=True)
sol=sp.solve(sp.Eq(Mk1*W,(Mk1-m)*Wp+m*(Wp-w)),Wp)[0]
print('increment',sp.simplify(sol-W))
assert sp.simplify(sol-W-m*w/Mk1)==0
# 一度に N 個
N=sp.symbols('N',positive=True,integer=True)
mc=sp.symbols('m_c',positive=True)
Vall=sp.solve(sp.Eq(0,(M-N*m)*W+N*m*(W-w)),W)[0]; print('all at once',Vall)
# 解答図の例：台車 m, 球 9 個 → M=10m
Ws=[Fr(0)]
for k in range(1,10):
    Ws.append(Ws[-1]+Fr(1,10-(k-1)))
for k in range(1,10):
    print(k,float(Ws[k]),' ball ground vel/w',float(Ws[k]-1))
print('W9=',float(Ws[9]),' all',9/10,' ln10',np.log(10))
first_forward=min(k for k in range(1,10) if Ws[k]>1); print('first forward ball',first_forward)
assert first_forward==7
# 不等式の確認（ランダム）
rng=np.random.default_rng(0)
for _ in range(1000):
    Mv=rng.uniform(2,50); mv=rng.uniform(0.01,1); n=int(rng.integers(2,int((Mv-0.5)/mv)+2)) if (Mv-0.5)/mv>2 else 2
    if n*mv>=Mv: continue
    s=sum(mv/(Mv-j*mv) for j in range(n))
    assert n*mv/Mv < s < np.log(Mv/(Mv-n*mv))
# 連続の極限
mu=sp.symbols('mu',positive=True); x=sp.symbols('x')
lim=sp.integrate(1/(M-x),(x,0,mu)); print('continuous',sp.simplify(lim))

# ---- (5) energy bookkeeping for the k-th launch ----
Wb=sp.symbols('W_b',nonnegative=True)   # W_{k-1}
Mp=Mk1
Wa=Wb+m*w/Mp
Ek=sp.simplify((Mp-m)*Wa**2/2+m*(Wa-w)**2/2-Mp*Wb**2/2)
print('(5a) E_k =',sp.factor(Ek)); assert sp.simplify(Ek-m*(Mp-m)*w**2/(2*Mp))==0
dK=(Mp-m)*(Wa**2-Wb**2)/2
ratio=sp.simplify(dK/Ek); print('(5b) dK/E =',ratio)
assert sp.simplify(ratio-(2*Wb/w+m/Mp))==0
# energy conservation: E = dK + dK_ball
dKball=m*(Wa-w)**2/2-m*Wb**2/2
assert sp.simplify(Ek-dK-dKball)==0
# ratio>1 <=> ball KE decreases
import random
for _ in range(2000):
    vals={m:random.uniform(.1,1),Mp:random.uniform(1.2,10),w:1.0,Wb:random.uniform(0,3)}
    r_=float(ratio.subs(vals)); db=float(dKball.subs(vals))
    assert (r_>1)==(db<0)
print('(5b) ratio>1 exactly when the ejected ball loses KE')
# example M=10m: ratio for each k
for k in range(1,10):
    Mpk=10-(k-1); print('   k=%d ratio=%.3f'%(k,2*float(Ws[k-1])+1/Mpk))
print('ok')
