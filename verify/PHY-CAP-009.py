# P58-1 共通の上極板をもつ2つのコンデンサー（直列）と誘電体板
import sympy as sp, numpy as np, math
from scipy.integrate import solve_ivp
e0,w,l,d,er,x,V,m=sp.symbols('epsilon0 w l d epsilon_r x V m',positive=True)
C0=e0*w*l/d   # 長さ l の空のコンデンサー1個の容量
CA=e0*w/d*((l/2+x)+er*(l/2-x))
CB=e0*w/d*((l/2-x)+er*(l/2+x))
# (1) 並列
print('(1) d(CA+CB)/dx =',sp.simplify(sp.diff(CA+CB,x)))
Cs=sp.simplify(CA*CB/(CA+CB))
Cs3=sp.simplify(Cs.subs(er,3))
print('(2) C(x) eps_r=3:',sp.factor(Cs3), ' / C0 =',sp.simplify(Cs3/C0))
Ueff=-sp.Rational(1,2)*Cs3*V**2
kV=sp.simplify(sp.diff(Ueff,x,2)); print('   k_V =',kV,' (C0V^2/l^2 =',sp.simplify(C0*V**2/l**2),')')
print('   d^3Ueff/dx^3=',sp.simplify(sp.diff(Ueff,x,3)))
Q=Cs3.subs(x,0)*V
UQ=Q**2/(2*Cs3)
kQ=sp.simplify(sp.diff(UQ,x,2).subs(x,0)); print('(3) k_Q at 0 =',kQ)
# 一般の er でも一致するか
CsG=Cs; QG=CsG.subs(x,0)*V
print('   general er: kV-kQ =',sp.simplify(sp.diff(-CsG*V**2/2,x,2).subs(x,0)-sp.diff(QG**2/(2*CsG),x,2).subs(x,0)))
# (4) x0 = l/2 から
vV=sp.sqrt(2*(Ueff.subs(x,l/2)-Ueff.subs(x,0))/m)
vQ=sp.sqrt(2*(UQ.subs(x,l/2)-UQ.subs(x,0))/m)
print('(4) vV=',sp.simplify(vV/sp.sqrt(C0/m)/V),'V sqrt(C0/m); vQ=',sp.simplify(vQ/sp.sqrt(C0/m)/V),' ratio',sp.simplify(vQ/vV))
# 電池の仕事 l/2 -> 0
Wb=V*(Cs3.subs(x,0)-Cs3.subs(x,l/2))*V
dUE=sp.Rational(1,2)*V**2*(Cs3.subs(x,0)-Cs3.subs(x,l/2))
print('   battery work',sp.simplify(Wb/C0/V**2),' dUE',sp.simplify(dUE/C0/V**2))
# 数値：電池接続時は厳密に単振動 → 周期 2π l sqrt(m/C0)/V
c0,ll,VV,mm=1.0,1.0,1.0,1.0
F=lambda xx: -c0*VV**2*xx/ll**2   # -dUeff/dx
so=solve_ivp(lambda t,y:[y[1],F(y[0])/mm],[0,20],[0.5,0],rtol=1e-12,atol=1e-12,dense_output=True)
t=np.linspace(0,20,200001); X=so.sol(t)[0]; zc=t[1:][(X[:-1]>0)&(X[1:]<=0)]
print('numeric period V',(zc[1]-zc[0]) if len(zc)>1 else None, ' theory',2*math.pi)
# 力を局所的に（端での引き込み力）求めて比較：A 内に入る向きの力 = 1/2 VA^2 (er-1) e0 w/d
xx=sp.symbols('xx')
VA=V*CB/(CA+CB); VB=V*CA/(CA+CB)
Floc=sp.Rational(1,2)*(er-1)*e0*w/d*(VB**2-VA**2)  # +x 向き（Bへ入る向き）
print('local force - (-dUeff/dx):',sp.simplify(Floc-(-sp.diff(-Cs*V**2/2,x))))
# 電荷一定の場合の数値積分（C0=V=l=m=1）：U = 1/(2(1-x^2))
FQ=lambda xx: -xx/(1-xx**2)**2
ev=lambda t,y:y[0]; ev.terminal=True; ev.direction=-1
so=solve_ivp(lambda t,y:[y[1],FQ(y[0])],[0,20],[0.5,0],events=ev,rtol=1e-12,atol=1e-12)
print('numeric vQ at center',abs(so.y_events[0][0][1]),' theory',1/math.sqrt(3))
