# M74-1 検証：凸レンズ＋液中の平面鏡（往復）系。光線行列（近軸）で独立計算
import sympy as sp
f,a,n,d,H=sp.symbols('f a n d H',positive=True)
def Tm(t): return sp.Matrix([[1,t],[0,1]])           # 空気換算距離 t の伝搬
def Lm(F): return sp.Matrix([[1,0],[-1/F,1]])
# 展開系：S -a- L -(H-d)- 液 d (換算 d/n) - 鏡 - 液 d - (H-d) - L - b -> 像
def system(av,nv,dv,Hv):
    gap=2*(Hv-dv)+2*dv/nv
    return Lm(f)*Tm(gap)*Lm(f)*Tm(av)
b=sp.symbols('b')
def image_dist(av,nv,dv,Hv):
    M=system(av,nv,dv,Hv)
    # 像点: (Tm(b)*M)[0,1]==0
    return sp.solve(sp.Eq((Tm(b)*M)[0,1],0),b)[0], M
# (1) 液体なし H=3f，像が S に一致：b=a
bb,_=image_dist(a,1,0,3*f)
print('(1) solutions a:',sp.solve(sp.Eq(bb,a),a))
# 液中で鏡の実効距離 D=H-d(1-1/n)
# (2) d=2f, H=3f，第2の一致位置 a'=5f/3 となる n
bb2,_=image_dist(a,n,2*f,3*f)
sol_all=sp.solve(sp.Eq(bb2,a),a); print('(2) coincidence a with liquid:',[sp.simplify(s) for s in sol_all])
print('    n from a=5f/3:',sp.solve(sp.Eq(bb2.subs(a,sp.Rational(5,3)*f),sp.Rational(5,3)*f),n))
# (3) n=4/3, a=9f/8
bb3,M=image_dist(sp.Rational(9,8)*f,sp.Rational(4,3),2*f,3*f)
print('(3) b =',sp.simplify(bb3))
Mt=sp.simplify(Tm(bb3)*M); print('    magnification A =',Mt[0,0])
# 中間像
b1=sp.Rational(9,8)*f*f/(sp.Rational(9,8)*f-f); print('    first image below lens:',b1,' 2D =',2*(3*f-2*f+2*f*sp.Rational(3,4)))
# 屈折を直接（スネル近軸ではなく実光線）で数値確認：小角度の実光線追跡
import numpy as np
fn=1.0; nn=4/3; dn=2.0; Hn=3.0; an=9/8
def trace(theta):
    # 展開して考えず，実光線を追う。y: 光軸からの距離，s: 進行方向に対する dy/dz
    y=an*np.tan(theta); s=np.tan(theta)      # S から下向きに出てレンズへ
    s=s-y/fn                                 # 薄レンズ
    y=y+s*(Hn-dn)                            # 空気中を液面まで
    s=np.tan(np.arcsin(np.sin(np.arctan(s))/nn))  # 液面で屈折（スネル則そのもの）
    y=y+s*dn                                 # 液中を鏡まで
    # 鏡で反射：進行方向が上向きに反転し，y の増え方（進行方向あたり）は不変
    y=y+s*dn                                 # 液中を液面まで
    s=np.tan(np.arcsin(np.sin(np.arctan(s))*nn))  # 液面で屈折して空気へ
    y=y+s*(Hn-dn)                            # レンズまで
    s=s-y/fn                                 # 薄レンズ
    return -y/s                              # 光軸と交わる点のレンズからの距離（上方）
print('    ray-trace b (small angles):',[trace(t) for t in [1e-4,1e-3,5e-3]])
# 一般の D での一致条件の因数分解
D=sp.symbols('D',positive=True)
b1=a*f/(a-f); a2=2*D-b1; b2=a2*f/(a2-f)
print('factor(b2-a):',sp.factor(sp.together(b2-a)))
