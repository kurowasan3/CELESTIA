# J126-1 検証：動かされる導体パイプ中を落下する磁石（相対速度で決まる電磁ブレーキとエネルギー収支）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,k,u,a,t=sp.symbols('m M g k u a t',positive=True)
v1=m*g/k
print('(1) v1=',v1,' P1=',sp.simplify(k*v1**2))
# (2) パイプ上向き速さ u：相対速さ w は m*0 = mg - k w -> w=v1。磁石の対地速度（下向き） = w-u
w=v1; Vdown=w-u
F=M*g+k*w           # パイプのつり合い：F - Mg - k w = 0
print('(2) Vdown=',Vdown,' F=',sp.simplify(F),' Phand=',sp.simplify(F*u),' Joule=',sp.simplify(k*w**2))
# エネルギー収支：手の仕事率 = パイプ位置E増加率 Mg u + 磁石位置E変化率 -mg Vdown + ジュール熱
print('  energy balance:',sp.simplify(F*u-(M*g*u-m*g*Vdown+k*w**2)))
# (3) u=v1 で静止。手の仕事率 - パイプ位置E増加率 = ジュール熱
print('(3) ',sp.simplify((F*u-M*g*u-k*w**2).subs(u,v1)))
# (4) 加速度 a で上昇：相対速度 wr（磁石がパイプに対して下向き）
wr=sp.Function('wr')
sol=sp.dsolve(sp.Eq(m*wr(t).diff(t),m*(g+a)-k*wr(t)),wr(t),ics={wr(0):m*g/k}).rhs
print('(4) w(t)=',sp.simplify(sol))
# 確認：対地座標で直接 磁石 m y'' = -mg + k(vp - y')（上向き正, vp=a t）
Fh=M*(g+a)+k*sol
print('    F(t)=',sp.simplify(Fh),' F(0)=',sp.simplify(Fh.subs(t,0)),' F(inf)=',sp.limit(Fh,t,sp.oo))
V=a*t-sol
print('    V(t)=',sp.expand(V))
# 数値積分（対地系）で照合
mn,gn,kn,an=0.2,9.8,4.0,3.0
def f(tt,y):  # y=[z_magnet, v_magnet] 上向き正
    vp=an*tt
    return [y[1], -gn + kn/mn*(vp-y[1])]
s=solve_ivp(f,[0,5],[0,-mn*gn/kn],rtol=1e-10,atol=1e-12,dense_output=True)
for tt in [0.5,1,3,5]:
    print(tt, s.sol(tt)[1], float(V.subs({a:an,g:gn,k:kn,m:mn,t:tt})))
# 長時間：V ≈ a(t - m/k) - mg/k
print('asym:',sp.simplify(sp.limit(V-(a*(t-m/k)-m*g/k),t,sp.oo)))
