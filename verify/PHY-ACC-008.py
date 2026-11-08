# J136-1 検証：理想変圧器による抵抗の変換と，最大電力（整合）と送電効率の対比
import sympy as sp, numpy as np
E,r,RL,a,k,Rl=sp.symbols('E r R_L a k R_l',positive=True)
# (1)(2) 発電機(E,r) - 変圧器 n1:n2 (a=n1/n2) - 負荷 RL
I1=E/(r+a**2*RL); P=sp.simplify(I1**2*a**2*RL)
print('(2) P(a)=',P)
sa=sp.solve(sp.diff(P,a),a); print('    a*=',sa)
Pm=sp.simplify(P.subs(a,sa[0])); print('    Pmax=',Pm,' eta=',sp.simplify(Pm/(E*I1.subs(a,sa[0]))))
# 別法：2次側から見たテブナン等価（起電力 E/a, 内部抵抗 r/a^2）
P2=(E/a)**2*RL/(r/a**2+RL)**2; print('    check thevenin:',sp.simplify(P2-P))
# (3) 昇圧1:k, 線路Rl, 降圧k:1, 負荷RL, 発電機内部抵抗0
Iline=k*E/(Rl+k**2*RL); PL=sp.simplify(Iline**2*k**2*RL); Pg=sp.simplify(E*k*Iline)
eta=sp.simplify(PL/Pg); print('(3) PL=',PL,' eta=',eta,' lim PL k->oo:',sp.limit(PL,k,sp.oo))
# 直接計算：負荷電流 IL = k Iline, 負荷電圧 = IL RL；発電機電流 = k Iline
# 線路損失 Iline^2 Rl と PL の和 = Pg
print('    balance:',sp.simplify(Pg-PL-Iline**2*Rl))
# (4) 数値 E=100, RL=1, Rl=2, eta>=0.98
kk=sp.symbols('kk',positive=True)
print('(4) k_min=',sp.solve(sp.Eq(eta.subs({RL:1,Rl:2,k:kk}),sp.Rational(98,100)),kk), ' sqrt98=',np.sqrt(98))
for kv in [9,10]:
    print('   k=',kv,' eta=',float(eta.subs({RL:1,Rl:2,k:kv})),' PL=',float(PL.subs({RL:1,Rl:2,k:kv,E:100})))
print('   k=1: eta=',float(eta.subs({RL:1,Rl:2,k:1})),' PL=',float(PL.subs({RL:1,Rl:2,k:1,E:100})))
# (5) 線路側から見た発電機・負荷：最大電力となる k（負荷の電力 PL(k) は k について単調増加か？）
print('   dPL/dk =',sp.factor(sp.diff(PL,k)))
