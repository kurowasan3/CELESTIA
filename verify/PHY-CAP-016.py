# T20242-1 エレクトレット発電：等価電池と操作の違いによる発熱
import sympy as sp
Q,q,z,h,d,eps,e0,S=sp.symbols('Q q z h d epsilon epsilon0 S',positive=True)
# 電場（上向き正）：誘電体中 (Q-q)/(eps S)，すき間 -q/(e0 S)
phi_top_minus_bottom = -((Q-q)*d/(eps*S) - q*z/(e0*S))
qz=sp.solve(sp.Eq(phi_top_minus_bottom,0),q)[0]
print('q(z)=',sp.simplify(qz))
E0=Q*d/(eps*S); Cd=eps*S/d; Cg=e0*S/z
print('phi = q(1/Cd+1/Cg)-E0 ?',sp.simplify(phi_top_minus_bottom-(q*(1/Cd+1/Cg)-E0)))
# エネルギー（電場のエネルギー）
U=lambda qq,zz: (Q-qq)**2*d/(2*eps*S)+qq**2*zz/(2*e0*S)
F=lambda qq: qq**2/(2*e0*S)   # 上電極に働く引力
# 電荷一定での引力 = -dU/dz
print('F=-dU/dz:',sp.simplify(-sp.diff(U(q,z),z)+F(q)))
qh=qz.subs(z,h)
# 操作P：z=0で閉(q=Q)→開いてhへ→閉→開いて0へ→閉
W_P=F(Q)*h - F(qh)*h
H1=U(Q,h)-U(qh,h) + 0   # 閉じたとき：電池仕事なし（外部回路は抵抗だけ）、q変化分のエネルギーが熱
H2=U(qh,0)-U(Q,0)
print('P: heat=work?',sp.simplify(H1+H2-W_P))
HP=sp.factor(sp.simplify(W_P))
Ch=1/(1/Cd+1/Cg.subs(z,h))
print('P heat formula check',sp.simplify(HP-(Q-qh)**2/2*(1/Ch+1/Cd)))
print('HP=',HP)
# 操作R：閉じたままゆっくりhへ→開いて0へ→閉
W_up=sp.integrate(F(qz),(z,0,h))
W_R=W_up - F(qh)*h
H_R=U(qh,0)-U(Q,0)
print('R: work=heat?',sp.simplify(W_R-H_R),'  H_R=(Q-qh)^2/(2Cd)?',sp.simplify(H_R-(Q-qh)**2/(2*Cd)))
# ゆっくり上げる間の熱が0であること：準静的ならエネルギー収支 W_up = U(qh,h)-U(Q,0)
print('quasi-static up: W_up = dU ?',sp.simplify(W_up-(U(qh,h)-U(Q,0))))
print('ratio HP/HR =',sp.simplify(HP/H_R))
# 数値例
import random
vals={Q:1.3,d:0.7,eps:2.1,e0:0.9,S:1.7,h:0.8}
print({k:float(v.subs(vals)) for k,v in [('HP',HP),('HR',H_R),('WP',W_P)]})
