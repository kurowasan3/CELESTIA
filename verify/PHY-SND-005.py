# M72-1 検証：真上通過の瞬間に聞く振動数（遅延時間），v と h の決定
import sympy as sp, numpy as np
from scipy.optimize import brentq
V,v,h,f0=sp.symbols('V v h f0',positive=True)
te=sp.symbols('t_e',real=True)
tr=te+sp.sqrt(v**2*te**2+h**2)/V          # 受信時刻（真上通過を t=0）
f=f0/sp.diff(tr,te)                        # 受信振動数
# (1) tr=0 となる te
sol=sp.solve(sp.Eq(tr,0),te); print('te at tr=0:',sol)
fstar=sp.simplify(f.subs(te,sol[0])); print('f* =',fstar)
print('x_e =',sp.simplify(v*sol[0]))
fmax=V/(V-v)*f0
r=sp.simplify(fstar/fmax); print('f*/fmax=',r)
vs=sp.solve(sp.Eq(r,sp.Rational(5,8)),v); print('v=',vs)
vv=sp.Rational(3,5)*V
fs=fstar.subs(v,vv)
# (3) 振動数 f*×8/17 を聞く時刻
target=fs*sp.Rational(8,17)
te2=sp.solve(sp.Eq(f.subs(v,vv),target),te); print('te for 8/17:',te2)
for t2 in te2:
    print(' tr=',sp.simplify(tr.subs({v:vv,te:t2})))
Vn=340.0; tau=5.0
hn=2*Vn*tau/5; print('h=',hn,'v=',0.6*Vn,'x at tau=',0.6*Vn*tau, '3h/2=',1.5*hn)
# 数値：独立に受信時刻を数値で求めて振動数を差分
vn=0.6*Vn; f0n=100.0
trf=lambda t: t+np.sqrt(vn**2*t**2+hn**2)/Vn
def f_at_reception(T):
    t=brentq(lambda x: trf(x)-T,-1e3,1e3); dt=1e-6
    return f0n*2*dt/(trf(t+dt)-trf(t-dt))
fs_n=f_at_reception(0.0); print('num f*/f0=',fs_n/f0n,'(25/16=',25/16,')')
print('num f(tau)/f* =',f_at_reception(tau)/fs_n,'(8/17=',8/17,')')
print('num f*/fmax =',fs_n/(f0n*Vn/(Vn-vn)))
