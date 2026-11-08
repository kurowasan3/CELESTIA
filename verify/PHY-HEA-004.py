# Q145-1: 放射シールドで囲んだ探査機の表面温度（片面だけ磨いたシールドの向きの効果）
import sympy as sp, numpy as np
I,A,sig,eps,N=sp.symbols('I A sigma epsilon N',positive=True)
Tp4,Ts4=sp.symbols('Tp4 Ts4',positive=True)
T04=I/(sig*A)
# (2) N 枚の黒いシールド: 各すき間で I = A sig (T_k^4 - T_{k+1}^4), 最外 I = A sig T_N^4
def black(Nn):
    Ts=sp.symbols('t0:%d'%(Nn+1))
    eqs=[sp.Eq(I,A*sig*(Ts[k]-Ts[k+1])) for k in range(Nn)]+[sp.Eq(I,A*sig*Ts[Nn])]
    return sp.solve(eqs,Ts,dict=True)[0][Ts[0]]
for Nn in (1,2,3,5): assert sp.simplify(black(Nn)-(Nn+1)*T04)==0
print('(2) T^4 = (N+1) T0^4  確認')
# (3) 内面黒・外面 ε
s3=sp.solve([sp.Eq(I,A*sig*(Tp4-Ts4)),sp.Eq(I,A*eps*sig*Ts4)],[Tp4,Ts4],dict=True)[0]
print('(3)',sp.simplify(s3[Tp4]/T04),sp.simplify(s3[Ts4]/T04))
# (4) 内面 ε・外面黒: 黒い探査機と灰色面の正味のやりとりを多重反射の級数で独立に導く
# 探査機が出す σTp^4 → 灰色面が ε 吸収，(1-ε) 反射 → 探査機（黒）がすべて吸収 → 終わり
k=sp.Symbol('k',integer=True,nonnegative=True)
absorbed_by_shield=eps*sig*Tp4               # 探査機の放射のうちシールドが吸収
emitted_by_shield_in=eps*sig*Ts4             # シールド内面の放射（探査機がすべて吸収）
net=absorbed_by_shield-emitted_by_shield_in
s4=sp.solve([sp.Eq(I,A*net),sp.Eq(I,A*sig*Ts4)],[Tp4,Ts4],dict=True)[0]
print('(4)',sp.simplify(s4[Tp4]/T04),sp.simplify(s4[Ts4]/T04))
assert sp.simplify(s3[Tp4]-s4[Tp4])==0
# 一般の2灰色面の式 (1/e1+1/e2-1)^-1 を多重反射の和で確かめ，e1=1 で上式に一致すること
e1,e2=sp.symbols('e1 e2',positive=True)
# 面1(e1)から出た放射のうち面2に吸収される割合：e2 Σ ((1-e1)(1-e2))^k
frac=e2/(1-(1-e1)*(1-e2))
q12=sp.simplify(e1*frac)   # σT1^4 あたり
print('2灰色面の有効係数',sp.simplify(q12-1/(1/e1+1/e2-1)))
assert sp.simplify(q12-1/(1/e1+1/e2-1))==0
assert sp.simplify(q12.subs(e1,1)-e2)==0
# 数値例（ε=0.05, T0=300 K）
T0=300.0; e=0.05
print('数値例 Tp=',T0*(1+1/e)**0.25,' Ts(3)=',T0*(1/e)**0.25,' Ts(4)=',T0, ' 黒1枚 Tp=',T0*2**0.25)
# 時間発展で定常状態に落ち着くことの確認（熱容量を適当に与えて数値積分）
from scipy.integrate import solve_ivp
sv=5.67e-8; Av=1.0; Iv=sv*Av*T0**4
for case in (3,4):
    def rhs(t,Y):
        Tp,Ts=Y
        if case==3: qin=sv*Av*(Tp**4-Ts**4); qout=e*sv*Av*Ts**4
        else: qin=e*sv*Av*(Tp**4-Ts**4); qout=sv*Av*Ts**4
        return [(Iv-qin)/1e4,(qin-qout)/1e3]
    so=solve_ivp(rhs,[0,5e6],[300,100],method='LSODA',rtol=1e-10)
    print('case',case,'定常 Tp,Ts =',so.y[0][-1],so.y[1][-1])
    assert abs(so.y[0][-1]-T0*(1+1/e)**0.25)<0.01
print('OK')
# ---- 新しい (5)：両面とも放射率 ε のシールド
# すき間（黒体の探査機と ε の面）の有効係数は q12(e1=1,e2=ε)=ε，外面は ε
s5=sp.solve([sp.Eq(I,A*q12.subs({e1:1,e2:eps})*sig*(Tp4-Ts4)),sp.Eq(I,A*eps*sig*Ts4)],[Tp4,Ts4],dict=True)[0]
print('(5) Tp^4/T0^4=',sp.simplify(s5[Tp4]/T04),' Ts^4/T0^4=',sp.simplify(s5[Ts4]/T04))
assert sp.simplify(s5[Tp4]/T04-2/eps)==0 and sp.simplify(s5[Ts4]/T04-1/eps)==0
ratio=sp.simplify(s5[Tp4]/s3[Tp4]); print('比 (5)/(3)=',ratio)
assert sp.simplify(ratio-2/(1+eps))==0 and sp.limit(ratio,eps,0)==2
for ev_ in (0.01,0.3,0.9): assert float(ratio.subs(eps,ev_))>1
# 時間発展でも同じ定常状態に落ち着く
def rhs5(t,Y):
    Tp,Ts=Y; qin=e*sv*Av*(Tp**4-Ts**4); qout=e*sv*Av*Ts**4
    return [(Iv-qin)/1e4,(qin-qout)/1e3]
so=solve_ivp(rhs5,[0,5e6],[300,100],method='LSODA',rtol=1e-10)
print('case5 定常 Tp,Ts =',so.y[0][-1],so.y[1][-1],' 理論',T0*(2/e)**0.25,T0*(1/e)**0.25)
assert abs(so.y[0][-1]-T0*(2/e)**0.25)<0.01 and abs(so.y[1][-1]-T0*(1/e)**0.25)<0.01
print('OK (5)')
