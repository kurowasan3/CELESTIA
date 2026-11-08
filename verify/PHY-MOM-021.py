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
# 数値例：台車 m, 球 9 個 → M=10m
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
print('ok')
