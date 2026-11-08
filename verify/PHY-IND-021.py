# T20122-1 : 斜めの境界をもつ磁場に進入するコンデンサー付き正方形回路
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp, dblquad
a,B,C,m,v,v0,X=sp.symbols('a B C m v v0 X',positive=True)
# 磁場領域 y<x と回路 [X-a,X]x[0,a] の重なり面積を数値積分で確認し，解析式と比較
def area_num(Xv,av=1.0):
    xs=np.linspace(Xv-av,Xv,20001); ys=np.clip(xs,0,av); return np.trapezoid(ys,xs)
def area_th(Xv,av=1.0):
    if Xv<=0: return 0
    if Xv<=av: return Xv**2/2
    if Xv<=2*av: return (av**2-(Xv-av)**2)/2+av*(Xv-av)
    return av**2
for Xv in [0.3,0.8,1.2,1.7,2.5]:
    print('X=%.1f area num %.6f th %.6f'%(Xv,area_num(Xv),area_th(Xv)))
A2=(a**2-(X-a)**2)/2+a*(X-a)
print('dA/dX (a<X<2a)=',sp.expand(sp.diff(A2,X)))
# (1)(2) 一定速度: Q=C*B*v*L, I=dQ/dt
L1=X; L2=2*a-X
I1=sp.diff(C*B*v*L1,X)*v; I2=sp.diff(C*B*v*L2,X)*v
print('I1=',I1,' I2=',I2)
# 外力 F_ext: 仕事率 F v = I * EMF
print('F1=',sp.simplify(I1*B*v*L1/v),' F2=',sp.simplify(I2*B*v*L2/v))
# 辺ごとのローレンツ力の直接計算 (反時計回り電流 I, B は +z)
def force_edges(Xv,Iv,av=1.0,Bv=1.0,N=4001):
    # 頂点 P(X-a,0),Q(X,0),R(X,a),S(X-a,a) 反時計回り P->Q->R->S->P
    pts=[(Xv-av,0),(Xv,0),(Xv,av),(Xv-av,av),(Xv-av,0)]
    F=np.zeros(2)
    for (x1,y1),(x2,y2) in zip(pts[:-1],pts[1:]):
        t=np.linspace(0,1,N); xs=x1+(x2-x1)*t; ys=y1+(y2-y1)*t
        inside=(ys<xs).astype(float)
        dl=np.array([x2-x1,y2-y1])/(N-1)
        # dF = I dl x B, B=(0,0,Bz): (dly*Bz, -dlx*Bz)
        F+=Iv*Bv*np.array([dl[1],-dl[0]])*inside.sum()
    return F
for Xv in [0.4,1.5]:
    print('X=',Xv,' F_mag per I B =',force_edges(Xv,1.0), ' expect x:',-(Xv if Xv<1 else 2-Xv))
# (3) 自由運動: エネルギー保存 (m+C B^2 L^2) v^2 = m v0^2 を数値積分と比較
mv,Cv,Bv,av,v0v=1.0,1.0,1.0,1.0,1.0   # C B^2 a^2 = m
def L(x): return x if 0<x<av else (2*av-x if av<=x<2*av else 0.0)
def Lp(x): return 1.0 if 0<x<av else (-1.0 if av<=x<2*av else 0.0)
def rhs(t,s):
    x,vv=s
    acc=-Cv*Bv**2*L(x)*Lp(x)*vv**2/(mv+Cv*Bv**2*L(x)**2)
    return [vv,acc]
sol=solve_ivp(rhs,[-0.5,6],[-0.5,v0v],max_step=1e-3,rtol=1e-10,atol=1e-12,dense_output=True)
for xq in [0.5,1.0,1.5,2.5]:
    k=np.argmin(abs(sol.y[0]-xq))
    print('X=%.2f v num %.6f  th %.6f'%(sol.y[0,k],sol.y[1,k],v0v*np.sqrt(mv/(mv+Cv*Bv**2*L(sol.y[0,k])**2))))
# 運動方程式（電荷 Q=CBvL を状態に入れない別の形）: m dv/dt = -I B L, I=dQ/dt を直接差分で
# 減速度最大の位置: |acc| = C B^2 m v0^2 X/(m+C B^2 X^2)^2
acc=C*B**2*m*v0**2*X/(m+C*B**2*X**2)**2
Xs=sp.solve(sp.diff(acc,X),X); print('X* =',Xs)
print('with CB^2a^2=m: X*=',sp.simplify(Xs[0].subs(C,m/(B**2*a**2))))
Qmax=C*B*v0*a*sp.sqrt(m/(m+C*B**2*a**2)); print('Qmax(CB^2a^2=m)=',sp.simplify(Qmax.subs(C,m/(B**2*a**2))))
print('acc max=',sp.simplify(acc.subs(C,m/(B**2*a**2)).subs(X,a/sp.sqrt(3))))
