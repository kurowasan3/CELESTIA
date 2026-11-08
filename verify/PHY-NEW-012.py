# J018-1 検証：水平に加速する水槽内の，底につないだ軽い球とふたからつるした重い球（見かけの重力と浮力の向き）
import sympy as sp, numpy as np
a,g,V,rho0,rho1,rho2,k=sp.symbols('a g V rho0 rho1 rho2 k',positive=True)
# 水槽の系：見かけの重力 g_eff = (-a, -g)（x は加速度の向き）
# 液中の圧力勾配 grad p = rho0 * g_eff  => 浮力 = -rho0 V g_eff = (rho0 V a, rho0 V g)
Tx,Ty=sp.symbols('Tx Ty')
# 軽い球（密度 rho1<rho0，底につなぐ）：浮力 + 重力 + 慣性力 + 張力 = 0
eqs=[rho0*V*a - rho1*V*a + Tx, rho0*V*g - rho1*V*g + Ty]
s1=sp.solve(eqs,[Tx,Ty]); print('light ball tension vector(on ball):',s1)
T1=sp.sqrt(s1[Tx]**2+s1[Ty]**2); print('|T1| =',sp.simplify(T1))
print('thread angle tan (from vertical), ball leans forward:',sp.simplify(-s1[Tx]/(-s1[Ty])))
# 重い球（密度 rho2>rho0，ふたからつるす）
eqs2=[rho0*V*a - rho2*V*a + Tx, rho0*V*g - rho2*V*g + Ty]
s2=sp.solve(eqs2,[Tx,Ty]); print('heavy ball tension vector:',s2, ' |T2|=',sp.simplify(sp.sqrt(s2[Tx]**2+s2[Ty]**2)))
# 液面の傾き tan = a/g
# 糸を切った軽い球の水槽に対する加速度（水の抵抗は無視）
acc=sp.simplify((rho0-rho1)/rho1*sp.sqrt(a**2+g**2)); print('relative accel after cut:',acc)
# (4) ばねの伸びの比 x_a/x_0 = sqrt(1+(a/g)^2) -> 5/4 なら a = 3g/4
print('a for ratio 5/4:',sp.solve(sp.Eq(sp.sqrt(1+(a/g)**2),sp.Rational(5,4)),a))
# 数値：見かけの重力の別解（地上系）: 軽い球は加速度 a で動く。地上系で F = m a
m1=rho1*V
eqx=sp.Eq(rho0*V*a + Tx, m1*a); eqy=sp.Eq(rho0*V*g - m1*g + Ty, 0)
print('ground frame check:',sp.solve([eqx,eqy],[Tx,Ty]))
# 数値で圧力場：水槽底の前端と後端の圧力差 rho0 a Lx
