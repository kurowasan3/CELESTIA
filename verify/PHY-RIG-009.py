# T20261-1 円周上を一定の接線加速度で加速する一輪車（棒の傾きとすべり出し）
import sympy as sp, numpy as np
m,g,a,R,t,mu=sp.symbols('m g a R t mu',positive=True)
v=a*t
# 慣性系：重心の加速度 = (接線 a, 中心向き v^2/R, 鉛直 0)。接触力 = m*acc + m g ez
Fc=sp.Matrix([m*a, m*v**2/R, m*g])   # (前, 中心向き, 上)
tan_theta=sp.sqrt(Fc[0]**2+Fc[1]**2)/Fc[2]
print('tan theta =',sp.simplify(tan_theta))
tan_psi=Fc[1]/Fc[0]; print('tan psi =',sp.simplify(tan_psi))
# 別解：棒とともに動く座標系で「見かけの重力」g_eff = g(下) + (-a 前) + (-v^2/R 中心向き) の逆向きに棒が沿う
geff=sp.Matrix([-a,-v**2/R,-g]); rod=-geff/geff.norm()
print('rod direction parallel to contact force:',sp.simplify(rod.cross(Fc)))
ts=sp.solve(sp.Eq(a**2+(a*t)**4/R**2,mu**2*g**2),t)
ts=[s for s in ts if s.is_positive!=False]
print('t_s candidates',ts)
tsv=sp.sqrt(R*sp.sqrt(mu**2*g**2-a**2))/a
print('check t_s:',sp.simplify((a**2+(a*tsv)**4/R**2)-mu**2*g**2))
print('tan theta at slip =',sp.simplify(tan_theta.subs(t,tsv)))
print('tan psi at slip =',sp.simplify(tan_psi.subs(t,tsv)))
s_=sp.simplify(a*tsv**2/2); print('distance s =',s_)
# 数値
vals={g:9.8,a:2.0,R:20.0,mu:0.5}
T=float(tsv.subs(vals)); print('t_s=',T,' theta(deg)=',np.degrees(np.arctan(float(tan_theta.subs(vals).subs(t,T)))),' atan(mu)=',np.degrees(np.arctan(0.5)))
# a=μg/√2 のとき psi=45度
print('psi at a=mu g/sqrt2:',sp.simplify(tan_psi.subs(t,tsv).subs(a,mu*g/sp.sqrt(2))))
