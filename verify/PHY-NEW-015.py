# Q080-1: 動ける三角台の粗い斜面上を滑る物体（コップの水面の傾き）
import sympy as sp, numpy as np
M,m,g,al,mu=sp.symbols('M m g alpha mu',positive=True)
A,N,ar=sp.symbols('A N a_r')  # A: 台の加速度(水平, +x), ar: 斜面に沿った相対加速度（下り向き正）
# 斜面は左上から右下へ下る：下り向き単位ベクトル t=(cos,-sin), 外向き法線 n=(sin,cos)
t=sp.Matrix([sp.cos(al),-sp.sin(al)]); n=sp.Matrix([sp.sin(al),sp.cos(al)])
def solve(sgn):
    # sgn=+1: 物体が下り向きに滑っている（摩擦は上り向き -t），sgn=-1: 上り向きに滑っている（摩擦は +t）
    fric=-sgn*mu*N*t
    R=N*n+fric                       # 台が物体に及ぼす力
    ablock=sp.Matrix([A,0])+ar*t
    eqs=list(m*ablock-(R+sp.Matrix([0,-m*g])))
    eqs.append(M*A+R[0])             # 台：水平方向（床はなめらか）
    s=sp.solve(eqs,[A,N,ar],dict=True)[0]
    Rs=R.subs(s); geff=-(Rs/m)       # 物体とともに動く系の見かけの重力 = g - a = -R/m
    # geff と斜面法線(-n)のなす角
    ang=sp.atan2(sp.simplify(geff.dot(t)),sp.simplify(-geff.dot(n)))
    return s,sp.simplify(ang)
sd,angd=solve(+1); su,angu=solve(-1)
print('down: A=',sp.simplify(sd[A]),' N=',sp.simplify(sd[N]),' a_r=',sp.simplify(sd[ar]))
print('down: geff tilt toward downslope = atan2(...)=',angd)
print('up  : A=',sp.simplify(su[A]),' a_r=',sp.simplify(su[ar]),' tilt=',angu)
# tilt = ±atan(mu) を数値で
vals={M:3,m:1,g:9.8,al:sp.pi/6,mu:sp.tan(sp.pi/18)}
print(float(angd.subs(vals)),float(angu.subs(vals)),float(sp.atan(vals[mu])))
assert abs(float(angd.subs(vals))-np.pi/18)<1e-12 and abs(float(angu.subs(vals))+np.pi/18)<1e-12
# M によらない：別の M でも
for Mv in [0.5,10,1000]:
    vv=dict(vals); vv[M]=Mv
    assert abs(float(angd.subs(vv))-np.pi/18)<1e-12
# 台の加速度の閉じた形
Ad=sp.simplify(sd[A]); target=m*g*sp.cos(al)*(sp.sin(al)-mu*sp.cos(al))/(M+m*sp.sin(al)*(sp.sin(al)-mu*sp.cos(al)))
print('A_down check:',sp.simplify(Ad-target), sp.simplify(Ad+target))
Nd=sp.simplify(sd[N]); print('N_down=',Nd)
# 水面と水平のなす角：上り α+λ，下り α−λ。β1=40°, β2=20° → α=30°, μ=tan10°
print('mu=',np.tan(np.radians(10)))
# 押す力がある場合：台に水平外力 F を加えても見かけの重力は -R/m のまま（Rの向きは法線から λ）
Fx=sp.symbols('F')
fric=-mu*N*t; R=N*n+fric; ablock=sp.Matrix([A,0])+ar*t
s=sp.solve(list(m*ablock-(R+sp.Matrix([0,-m*g])))+[M*A+R[0]-Fx],[A,N,ar],dict=True)[0]
geff=-(R.subs(s)/m); ang=sp.simplify(sp.atan2(sp.simplify(geff.dot(t)),sp.simplify(-geff.dot(n))))
print('with F: tilt =',ang)
print('OK')
